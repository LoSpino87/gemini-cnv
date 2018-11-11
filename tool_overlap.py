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
from compression import pack_blob

def _sample_name(database):
	names = []
	query = "SELECT name FROM samples"
	name = GeminiQuery.GeminiQuery(database)
	name.run(query)
	for n in name:
		names.append(str(n))
	return names

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
		sel_sample = args.sample.split(',')
		gts = 'gts.' + ', gts.'.join([str(s) for s in sel_sample])
		gt_filter  = " != './.' or".join([s for s in gts.split(',')]) + " != './.' "

		args.query = """select chrom, start, end, alt, variant_id,""" + gts + """from variants_cnv"""
	else :
		gt_filter = None
		sel_sample = _sample_name(database = args.db)
		args.query = """select chrom, start, end, alt, variant_id,(gts).(*) from variants_cnv"""

	VAR = GeminiQuery.GeminiQuery(args.db)
	VAR.run(args.query,gt_filter)
	var_string = ""
	for i in VAR:
		print i
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
		query_gen = "SELECT (gts).(*), (gt_types).(*),(gt_phases).(*),(gt_depths).(*),(gt_ref_depths).(*),(gt_alt_depths).(*),\
				(gt_quals).(*),(gt_copy_numbers).(*),(gt_phred_ll_homref).(*),(gt_phred_ll_het).(*),(gt_phred_ll_homalt).(*) \
				from variants_cnv v where v.variant_id=={}".format(str(row[4]))
		geno = GeminiQuery.GeminiQuery(args.db)
		geno.run(query_gen)
		for v in geno:
			over_res = dict(chrom_A=row[0],
							start_A=int(row[1]),
							end_A=int(row[2]),
							alt=row[3],
							gts=pack_blob(v['gts']),
							gt_types = pack_blob(v['gt_types']),
							gt_phases = pack_blob(v['gt_phases']) ,
		                    gt_depths=pack_blob(v['gt_depths']),
							gt_ref_depths=pack_blob(v['gt_ref_depths']),
		                    gt_alt_depths=pack_blob(v['gt_alt_depths']),
		                    gt_quals=pack_blob(v['gt_quals']),
							gt_copy_numbers=pack_blob(v['gt_copy_numbers']),
							gt_phred_ll_homref=pack_blob(v['gt_phred_ll_homref']),
							gt_phred_ll_het=pack_blob(v['gt_phred_ll_het']),
							gt_phred_ll_homalt=pack_blob(v['gt_phred_ll_homalt']),
							chrom_B = row[-15],
							start_B = int(row[-14]),
							end_B = int(row[-13]),
							type = row[-12],
							num_variants = row[-11],
							num_samples = row[-10],
							African = row[-9],
							Asian = row[-8],
							European = row[-7],
							Mexican = row[-6],
							Middle_east = row[-5],
							Native_american = row[-4],
							Oceania = row[-3],
							South_american= row[-2],
							overlap_bp = int(row[-1])+1,
							len_A = int(row[2])-int(row[1])+1,
							len_B = int(row[-13])-int(row[-14])+1,
							overlap_A_perc = str(round(float(int(row[-1])+1)/(int(row[2])-int(row[1])+1)*100,1)),
							overlap_B_perc = str(round(float(int(row[-1])+1)/(int(row[-13])-int(row[-14])+1)*100,1)))

			if over_res['overlap_A_perc'] == '100.0': over_res['overlap_A_perc'] = '100'
			if over_res['overlap_A_perc'] == '0.0': over_res['overlap_A_perc'] = '<0.1'
			if over_res['overlap_B_perc'] == '100.0': over_res['overlap_B_perc'] = '100'
			if over_res['overlap_B_perc'] == '0.0': over_res['overlap_B_perc'] = '<0.1'
			over_res['jaccard_index'] = jaccard_index(over_res['overlap_bp'],over_res['len_A'],over_res['len_B'])
			var_and_cnv_b.append(over_res)

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
		r['uid']=id
		overlap_result.append(r)

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
		query_gen = "SELECT (gts).(*), (gt_types).(*),(gt_phases).(*),(gt_depths).(*),(gt_ref_depths).(*),(gt_alt_depths).(*),\
				(gt_quals).(*),(gt_copy_numbers).(*),(gt_phred_ll_homref).(*),(gt_phred_ll_het).(*),(gt_phred_ll_homalt).(*) \
				from variants_cnv v where v.variant_id=={}".format(row[4])
		geno = GeminiQuery.GeminiQuery(args.db)
		geno.run(query_gen)
		for v in geno:
			over_res = dict(chrom_A=row[0],
							start_A=int(row[1]),
							end_A=int(row[2]),
							alt=row[3],
							gts=pack_blob(v['gts']),
							gt_types = pack_blob(v['gt_types']),
							gt_phases = pack_blob(v['gt_phases']) ,
		                    gt_depths=pack_blob(v['gt_depths']),
							gt_ref_depths=pack_blob(v['gt_ref_depths']),
		                    gt_alt_depths=pack_blob(v['gt_alt_depths']),
		                    gt_quals=pack_blob(v['gt_quals']),
							gt_copy_numbers=pack_blob(v['gt_copy_numbers']),
							gt_phred_ll_homref=pack_blob(v['gt_phred_ll_homref']),
							gt_phred_ll_het=pack_blob(v['gt_phred_ll_het']),
							gt_phred_ll_homalt=pack_blob(v['gt_phred_ll_homalt']),
							chrom_B = row[-5],
							start_B = int(row[-4]),
							end_B = int(row[-3]),
							opt_field = row[-2],
							overlap_bp = int(row[-1])+1,
							len_A = int(row[2])-int(row[1])+1,
							len_B = int(row[-3])-int(row[-4])+1,
							overlap_A_perc = str(round(float(int(row[-1])+1)/(int(row[2])-int(row[1])+1)*100,1)),
							overlap_B_perc = str(round(float(int(row[-1])+1)/(int(row[-3])-int(row[-4])+1)*100,1)))


			if over_res['overlap_A_perc'] == '100.0': over_res['overlap_A_perc'] = '100'
			if over_res['overlap_A_perc'] == '0.0': over_res['overlap_A_perc'] = '<0.1'
			if over_res['overlap_B_perc'] == '100.0': over_res['overlap_B_perc'] = '100'
			if over_res['overlap_B_perc'] == '0.0': over_res['overlap_B_perc'] = '<0.1'
			over_res['jaccard_index'] = jaccard_index(over_res['overlap_bp'],over_res['len_A'],over_res['len_B'])
			var_and_cnv_b.append(over_res)

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
		r['uid']=id
		overlap_custom_result.append(r)

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
		if str(i['alt_A']) == args.alt_par.upper():
			overlap_filt_a.append(i)
	return overlap_filt_a

def overlap_filt_len(args,result):
	overlap_filt_l = []
	for i in result:
		if args.int_len_min:
			if i['overlap_bp'] > int(args.int_len_min):
				if args.int_len_max:
					if i['overlap_bp'] < int(args.int_len_max):
						overlap_filt_l.append(i)
				else:
					overlap_filt_l.append(i)
		else:
			if args.int_len_max:
				if i['overlap_bp'] < int(args.int_len_max):
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

	var_no_cnv_b=[]
	id=0
	for row in var_no_cnv:
		query_gen = "SELECT (gts).(*), (gt_types).(*),(gt_phases).(*),(gt_depths).(*),(gt_ref_depths).(*),(gt_alt_depths).(*),\
				(gt_quals).(*),(gt_copy_numbers).(*),(gt_phred_ll_homref).(*),(gt_phred_ll_het).(*),(gt_phred_ll_homalt).(*) \
				from variants_cnv v where v.variant_id=={}".format(row[4])
		geno = GeminiQuery.GeminiQuery(args.db)
		geno.run(query_gen)
		for v in geno:
			over_res = dict(chrom=row[0],
							start=int(row[1]),
							end=int(row[2]),
							alt=row[3],
							gts=pack_blob(v['gts']),
							gt_types = pack_blob(v['gt_types']),
							gt_phases = pack_blob(v['gt_phases']) ,
		                    gt_depths=pack_blob(v['gt_depths']),
							gt_ref_depths=pack_blob(v['gt_ref_depths']),
		                    gt_alt_depths=pack_blob(v['gt_alt_depths']),
		                    gt_quals=pack_blob(v['gt_quals']),
							gt_copy_numbers=pack_blob(v['gt_copy_numbers']),
							gt_phred_ll_homref=pack_blob(v['gt_phred_ll_homref']),
							gt_phred_ll_het=pack_blob(v['gt_phred_ll_het']),
							gt_phred_ll_homalt=pack_blob(v['gt_phred_ll_homalt']))

			over_res['len'] = over_res['end']-over_res['start']+1
			var_no_cnv_b.append(over_res)

	# length filter
	if args.int_len_min or args.int_len_max:
		var_no_cnv_b = no_overlap_filt_len(args=args,result=var_no_cnv_b)
	# alteration filter
	if args.alt_par:
		var_no_cnv_b = overlap_filt_alt(args=args,result=var_no_cnv_b)

	overlap_result = []
	for r in var_no_cnv_b:
		id += 1
		r['uid']=id
		overlap_result.append(r)

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
			if i['len'] > int(args.int_len_min):
				if args.int_len_max:
					if i['len'] < int(args.int_len_max):
						overlap_filt_l.append(i)
				else:
					overlap_filt_l.append(i)
		else:
			if args.int_len_max:
				if i['len'] < int(args.int_len_max):
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
			print 'There is no map loaded. Please select a cnv map from DGV (--dgv_cnvmap) or a BED (--bed).'

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
