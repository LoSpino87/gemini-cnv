#A tool for valutate overlapping of genome
#!/usr/bin/env python

# Python imports
import os.path
import re
import sys
import pybedtools
from pybedtools import BedTool
import sqlalchemy as sql
from sqlalchemy.orm import create_session

# Gemini imports
import GeminiQuery
import database
import dgv_table
import export_vcf

# get variants from cnv map file
def get_dgv_map(args):
	"""
	Define the dgv_map table
	"""
	c, metadata = database.get_session_metadata(args.db)
	# drop table if already exists
	c.execute("DROP TABLE if exists dgv_map")
	# create table
	database.create_dgv_table(c,metadata,args)
	#unique identifier for each entry
	i = 0
	contents = dgv_map_c = []

	with open(args.dgv_cnvmap,'r') as g:
		next(g)
		for line in g:
			col = line.strip().split("\t")
			table = dgv_table.dgv_map(col)
			i += 1
			dgv_map_c = [str(i), table.chr, table.start, table.end, table.state, table.id, table.type, table.num_variants,
                  table.num_samples,table.num_samples_multicounted, table.num_studies, table.variants, table.samples,
                  table.studies, table.African, table.Asia, table.European, table.Mexican, table.Middle_East, table.Native_American,
                  table.Oceania, table.South_American]
			contents.append(dgv_map_c)
			if i % 1000 == 0:
				database.insert_dgv_map(c,metadata, contents)
				contents = []
		database.insert_dgv_map(c, metadata, contents)
		database.close_and_commit(c)

def get_cnv_map(args):
	"""
	Define and populate cnv_map table
	"""

	c, metadata = database.get_session_metadata(args.db)
	# drop table if already exists
	c.execute("DROP TABLE if exists cnv_map")
	#create table
	database.create_cnv_map(c,metadata,args)
	#unique identifier for each entry
	i = 0
	contents = cnv_map_c = []
	with open(args.bed_cnvmap,'r') as g:
		next(g)
		for line in g:
			col = line.strip().split("\t")
			table = dgv_table.cnv_custom_map(col)
			i += 1
			cnv_map_c = [str(i), table.chr, table.start, table.end, table.opt_field]
			contents.append(cnv_map_c)
			if i % 1000 == 0:
				database.insert_cnv_map(c,metadata, contents)
				contents = []
		database.insert_cnv_map(c, metadata, contents)
		database.close_and_commit(c)

def overlap_main(args):
	if args.dgv_cnvmap:
		overlap_res(args)
		print_rec(args)
	if args.bed_cnvmap:
		overlap_custom_res(args)
		print_custom_rec(args)

def extract_var(args):
	# extract data from variants table and create relative BED object
	if args.sample:
		gt_col = 'gts.' + ', gts.'.join([s for s in str(args.sample).split(',')])
		gt_filter  = " != './.' or".join([s for s in gt_col.split(',')]) + " != './.' "
		args.query = """select chrom, start, end, alt, """ + gt_col + """ from variants_cnv"""
	else :
		gt_filter = None
		args.query = 'select chrom, start, end, alt from variants_cnv'

	VAR = GeminiQuery.GeminiQuery(args.db)
	VAR.run(args.query,gt_filter)
	var_string = ""
	for i in VAR:
		if args.sample:
			i = "\t".join([x for x in str(i).split('\t')[:4]])

		var_string = var_string + "\n" + str(i)

	var_bed = pybedtools.BedTool(var_string, from_string=True)
	return var_bed

def extract_cnv(args):
	# extract data from dgv_map table and create relative pyBED object
	args.query = """select chr, start, end, type, num_variants, num_samples,
					African, Asian, European, Mexican, Middle_east,
					Native_american, Oceania, South_american
					from dgv_map"""
	CNV = GeminiQuery.GeminiQuery(args.db)
	CNV.run(args.query)
	cnv_string = ""
	for i in CNV:
		cnv_string = cnv_string + "\n" + str(i)
	cnv_bed = pybedtools.BedTool(cnv_string, from_string=True)
	return cnv_bed

def extract_cnv_custom(args):
	# extract data from dgv_map table and create relative pyBED object
	args.query = "select chr, start, end, opt_field FROM cnv_map"
	CNV = GeminiQuery.GeminiQuery(args.db)
	CNV.run(args.query)
	cnv_string = ""
	for i in CNV:
		cnv_string = cnv_string + "\n" + str(i)
	cnv_bed = pybedtools.BedTool(cnv_string, from_string=True)
	return cnv_bed

def overlap_res(args):
	# get cnv map
	get_dgv_map(args)

	# extract data
	var_bed = extract_var(args)
	cnv_bed = extract_cnv(args)

	if args.v:
		no_overlap(args,var_bed,cnv_bed)
		#if args.out_file:
		#	export_vcf.export_vcf('no_overlap', args.db, args.out_file)
	else:
		overlap(args,var_bed,cnv_bed)
		#if args.out_file:
		#	export_vcf.export_vcf('overlap', args.db, args.out_file)

def overlap_custom_res(args):
	# get cnv custom map
	get_cnv_map(args)

	# extract data
	var_bed = extract_var(args)
	cnv_bed = extract_cnv_custom(args)

	if args.v:
		no_overlap(args, var_bed, cnv_bed)
		#if args.out_file:
		#	export_vcf.export_vcf('no_overlap', args.db, args.out_file)
	else:
		overlap_custom(args,var_bed,cnv_bed)
		#if args.out_file:
		#	export_vcf.export_vcf('overlap', args.db, args.out_file)

def overlap(args,var_bed,cnv_bed):
	if args.f_par:
		if args.r:
			var_and_cnv = overlap_reciprocal(args,var_bed,cnv_bed)
		else:
			var_and_cnv = overlap_fraction(args,var_bed,cnv_bed)
	else:
		var_and_cnv = var_bed.intersect(cnv_bed, wo= True)

	#show result
	var_and_cnv_b = []
	for row in var_and_cnv:
		inter_overlap = int(row[18])+1
		length_A = int(row[2])-int(row[1])+1
		length_B = int(row[6])-int(row[5])+1
		perc_A = str(round(float(inter_overlap)/length_A*100,1))
		if perc_A == '100.0': perc_A = '100'
		if perc_A == '0.0': perc_A = '<0.1'
		perc_B = str(round(float(inter_overlap)/length_B*100,1))
		if perc_B == '100.0': perc_B = '100'
		if perc_B == '0.0': perc_B = '<0.1'
		j_index = jaccard_index(inter_overlap,length_A,length_B)

		var_and_cnv_b.append([row[0],row[1],row[2],length_A,perc_A,row[3],row[4],row[5],row[6],length_B,perc_B,row[7],inter_overlap,j_index,row[8],row[9],row[10],row[11],row[12],row[13],row[14],row[15],row[16],row[17]])

	# length filter
	if args.int_len_min or args.int_len_max:
		var_and_cnv_b = overlap_filt_len(args=args,result=var_and_cnv_b)
	# alteration filter
	if args.alt_par:
		var_and_cnv_b = overlap_filt_alt(args=args,result=var_and_cnv_b)

	# insert id for each row of table
	overlap_result = []
	id = 0
	for r in var_and_cnv_b:
		id += 1
		overlap_result.append([id,r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],r[10],r[11],r[12],r[13],r[14],r[15],r[16],r[17],r[18],r[19],r[20],r[21],r[22],r[23]])

	e = sql.create_engine(database.get_path(args.db), isolation_level=None)
	e.connect().connection.connection.text_factory = str
	metadata = sql.MetaData(bind=e)
	session = create_session(bind=e, autocommit=False, autoflush=False)
	over_table = sql.Table('overlap', metadata)

	c, metadata = database.get_session_metadata(args.db)
	c.execute("DROP TABLE if exists overlap")
	database.create_overlap_result(c,metadata,args)
	database.insert_overlap(c, metadata,overlap_result)
	database.close_and_commit(c)

def overlap_custom(args,var_bed,cnv_bed):
	if args.f_par:
		if args.r:
			var_and_cnv = overlap_reciprocal(args,var_bed,cnv_bed)
		else:
			var_and_cnv = overlap_fraction(args,var_bed,cnv_bed)
	else:
		var_and_cnv = var_bed.intersect(cnv_bed, wo= True)

	#show result
	var_and_cnv_b = []

	for row in var_and_cnv:
		inter_overlap = int(row[8])+1
		length_A = int(row[2])-int(row[1])+1
		length_B = int(row[6])-int(row[5])+1
		perc_A = str(round(float(inter_overlap)/length_A*100,1))
		if perc_A == '100.0': perc_A = '100'
		if perc_A == '0.0': perc_A = '<0.1'
		perc_B = str(round(float(inter_overlap)/length_B*100,1))
		if perc_B == '100.0': perc_B = '100'
		if perc_B == '0.0': perc_B = '<0.1'
		j_index = jaccard_index(inter_overlap,length_A,length_B)
		var_and_cnv_b.append([row[0],row[1],row[2],length_A,perc_A,row[3],row[4],row[5],row[6],length_B,perc_B,row[7],inter_overlap,j_index])

	# length filter
	if args.int_len_min or args.int_len_max:
		var_and_cnv_b = overlap_filt_len(args=args,result=var_and_cnv_b)
	# alteration filter
	if args.alt_par:
		var_and_cnv_b = overlap_filt_alt(args=args,result=var_and_cnv_b)

	# insert id for each row of table
	overlap_custom_result = []
	id = 0

	for r in var_and_cnv_b:
		id += 1
		overlap_custom_result.append([id,r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],r[10],r[11],r[12],r[13]])

	e = sql.create_engine(database.get_path(args.db), isolation_level=None)
	e.connect().connection.connection.text_factory = str
	metadata = sql.MetaData(bind=e)
	session = create_session(bind=e, autocommit=False, autoflush=False)
	over_table = sql.Table('overlap_custom', metadata)

	c, metadata = database.get_session_metadata(args.db)
	c.execute("DROP TABLE if exists overlap_custom")
	database.create_overlap_custom_result(c,metadata,args)
	database.insert_overlap_custom(c, metadata,overlap_custom_result)
	database.close_and_commit(c)

def overlap_fraction(args,var_bed,cnv_bed):
	var_and_cnv = var_bed.intersect(cnv_bed, f = args.f_par, wo = True)
	return var_and_cnv

def overlap_reciprocal(args,var_bed,cnv_bed):
	var_and_cnv = var_bed.intersect(cnv_bed, f = args.f_par, r = True, wo = True)
	return var_and_cnv

# filtering result
def overlap_filt_alt(args,result):
	overlap_filt_a = []
	for i in result:
		if str(i[5]) == args.alt_par.upper():
			overlap_filt_a.append(i)
	return overlap_filt_a

def overlap_filt_len(args,result):
	overlap_filt_l = []
	for i in result:
		if args.int_len_min:
			if i[12] > int(args.int_len_min):
				if args.int_len_max:
					if i[12] < int(args.int_len_max):
						overlap_filt_l.append(i)
				else:
					overlap_filt_l.append(i)
		else:
			if args.int_len_max:
				if i[12] < int(args.int_len_max):
					overlap_filt_l.append(i)
			else:
				overlap_filt_l.append(i)
	return overlap_filt_l

# add jaccard index
def jaccard_index(interval,seq1,seq2):
    index = round(float(float(interval)/(float(seq1) + float(seq2))),3)
    return index

# funtion for no overlap variants
def no_overlap(args,var_bed,cnv_bed):
	if args.f_par:
		var_no_cnv = var_bed.intersect(cnv_bed, f = args.f_par, v = True)
		if args.r:
			var_no_cnv = var_bed.intersect(cnv_bed, f = args.f_par, r= True, v = True)
	else:
		var_no_cnv = var_bed.intersect(cnv_bed, v = True)

	overlap_result = []
	id = 0
	for r in var_no_cnv:
		id += 1
		len = int(r[2])-int(r[1])+1
		overlap_result.append([id,r[0],r[1],r[2],len,r[3]])

	# length filter
	if args.int_len_min or args.int_len_max:
		overlap_result = no_overlap_filt_len(args=args,result=overlap_result)
	# alteration filter
	if args.alt_par:
		overlap_result = overlap_filt_alt(args=args,result=overlap_result)

	e = sql.create_engine(database.get_path(args.db), isolation_level=None)
	e.connect().connection.connection.text_factory = str
	metadata = sql.MetaData(bind=e)
	session = create_session(bind=e, autocommit=False, autoflush=False)
	over_table = sql.Table('no_overlap', metadata)

	c, metadata = database.get_session_metadata(args.db)
	c.execute('''DROP TABLE if exists no_overlap''')
	database.create_no_overlap_result(c,metadata,args)
	database.insert_no_overlap(c, metadata, overlap_result)
	database.close_and_commit(c)

def no_overlap_filt_len(args,result):
	overlap_filt_l = []
	for i in result:
		if args.int_len_min:
			if i[4] > int(args.int_len_min):
				if args.int_len_max:
					if i[4] < int(args.int_len_max):
						overlap_filt_l.append(i)
				else:
					overlap_filt_l.append(i)
		else:
			if args.int_len_max:
				if i[4] < int(args.int_len_max):
					overlap_filt_l.append(i)
			else:
				overlap_filt_l.append(i)
	return overlap_filt_l

# run tool from main parser
def run(parser, args):
	if os.path.exists(args.db):
		if args.dgv_cnvmap or args.bed_cnvmap:
			overlap_main(args)
		else:
			print 'There is no map laoded. Please select a cnv map from DGV (--dgv_cnvmap) or a BED (--bed).'

# name utility for browser
def name_dgv(database):
	name = "SELECT resource FROM resources WHERE name='dgv_cnvmap'"
	nm = GeminiQuery.GeminiQuery(database)
	nm.run(name)
	for n in nm:
		return n

# print result from DGV ovelap result table
def print_rec(args):
	if args.v == True:
		args.query = "select * from no_overlap"
	else:
		args.query = "select * from overlap"

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(args.query)
	print res.header
	for row in res:
		print row

# print result from CNV custom map overlap table
def print_custom_rec(args):
	if args.v == True:
		args.query = "select * from no_overlap"
	else:
		args.query = "select * from overlap_custom"

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(args.query)
	print res.header
	for row in res:
		print row
