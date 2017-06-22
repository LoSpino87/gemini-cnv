#A tool for valutate overlapping of genome
#!/usr/bin/env python

# Python imports
import os.path
import re
import sys
import pybedtools
from pybedtools import BedTool
import sqlalchemy as sql
from sqlalchemy.orm import mapper, create_session
import math

# Gemini imports
import GeminiQuery
import database


def run(parser, args):
	if os.path.exists(args.db):
		overlap_main(args)

def overlap_main(args):
	overlap(args)
	args.query = "select * from overlap"

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(args.query)
	for row in res:
		print row

def overlap(args):
	# extract data from variants table and create relative BED object
	args.query = 'select chrom, start, end, alt from variants_cnv'
	VAR = GeminiQuery.GeminiQuery(args.db)
	VAR.run(args.query)
	var_string = ""
	for i in VAR:
		var_string = var_string + "\n" + str(i)
	var_bed = pybedtools.BedTool(var_string, from_string=True)


	# extract data from dgv_map table and create relative BED object
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

	#overlap
	var_and_cnv = []
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
	metadata.reflect(bind=e)
	metadata.create_all()
	session = create_session(bind=e, autocommit=False, autoflush=False)

	database.empty_overlap_table(session=session,metadata=metadata)
	database.insert_overlap(session=session, metadata=metadata, overlap = overlap_result)

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

def jaccard_index(interval,seq1,seq2):
    index = round(float(float(interval)/(float(seq1) + float(seq2))),3)
    return index

def overlap_fraction(args,var_bed,cnv_bed):
	var_and_cnv = var_bed.intersect(cnv_bed, f = args.f_par, wo = True)
	return var_and_cnv

def overlap_reciprocal(args,var_bed,cnv_bed):
	var_and_cnv = var_bed.intersect(cnv_bed, f = args.f_par, r = True, wo = True)
	return var_and_cnv
