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

# Gemini imports
import GeminiQuery
import database


def run(parser, args):
	if os.path.exists(args.db):
		overlap_main(args)

def overlap_main(args):
    overlap(args)
    args.query = 'select * from overlap'
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
	args.query = 'select chr, start, end, type from dgv_map'
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
	id = 0
	for row in var_and_cnv:
		id += 1
		inter_overlap = int(row[8])+1
		length_A = int(row[2])-int(row[1])+1
		length_B = int(row[6])-int(row[5])+1
		perc_A = round(float(inter_overlap)/length_A*100,1)
		if perc_A == 100.0: perc_A = 100
		if perc_A == 0.0: perc_A = '<0.1'
		perc_B = round(float(inter_overlap)/length_B*100,1)
		if perc_B == 100.0: perc_B = 100
		if perc_B == 0.0: perc_B = '<0.1'
		j_index = jaccard_index(inter_overlap,length_A,length_B)
		var_and_cnv_b.append([id,row[0],row[1],row[2],length_A,perc_A,row[3],row[4],row[5],row[6],length_B,perc_B,row[7],inter_overlap,j_index])

	e = sql.create_engine(database.get_path(args.db), isolation_level=None)
	e.connect().connection.connection.text_factory = str
	metadata = sql.MetaData(bind=e)
	metadata.reflect(bind=e)
	metadata.create_all()
	session = create_session(bind=e, autocommit=False, autoflush=False)

	database.empty_overlap_table(session=session,metadata=metadata)
	database.insert_overlap(session=session, metadata=metadata, overlap = var_and_cnv_b)

def jaccard_index(interval,seq1,seq2):
    index = round(float(float(interval)/(float(seq1) + float(seq2))),3)
    return index

def overlap_fraction(args,var_bed,cnv_bed):
	var_and_cnv = var_bed.intersect(cnv_bed, f = args.f_par, wo = True)
	return var_and_cnv

def overlap_reciprocal(args,var_bed,cnv_bed):
	var_and_cnv = var_bed.intersect(cnv_bed, f = args.f_par, r = True, wo = True)
	return var_and_cnv
