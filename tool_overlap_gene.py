#A tool for valutate overlapping of genome
#!/usr/bin/env python

# Python imports
import os.path
import os
import re
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sqlalchemy as sql
import numpy as np
import math
import sqlalchemy as sql
from sqlalchemy.orm import create_session

# Gemini imports
import GeminiQuery
import gene_table
import database
from compression import pack_blob

gene = []
alt = np.array([])
samples = []

def run(parser, args):
	if os.path.exists(args.db):
		samples = _sample_name(database = args.db)

		if args.gene_map:
			overlap_custom_gene(args)
		else:
			overlap_gene_main(args)

		print_rec(args)

def _extract_data(alt,row,sel_sample,smp2idx):
	temp = []
	for s in sel_sample:
		idx = smp2idx[s]
		if row['gts'].take(idx) != './.':
			if row['alt'] == 'DUP':
				temp.append(1)
			elif row['alt'] == 'DEL':
				temp.append(-1)
			else:
				temp.append(0)
		else:
			temp.append(0)
	alt = np.concatenate((alt,[temp]),axis=0)
	return alt

def _genename(row):
	gene = str(row['gene'])
	if gene == 'None':
		gene = str(row['synonym'])
		if gene == 'None':
			gene = str(row['ensembl_gene_id'])
	return gene

def _query(args):

	if args.sample:
		sel_sample = args.sample.split(',')

		gts = 'gts.' + ', gts.'.join([str(s) for s in sel_sample])
		gt_types = 'gt_types.' + ', gt_types.'.join([str(s) for s in sel_sample])
		gt_phases = 'gt_phases.' + ', gt_phases.'.join([str(s) for s in sel_sample])
		gt_depths = 'gt_depths.' + ', gt_depths.'.join([str(s) for s in sel_sample])
		gt_ref_depths = 'gt_ref_depths.' + ', gt_ref_depths.'.join([str(s) for s in sel_sample])
		gt_alt_depths = 'gt_alt_depths.' + ', gt_alt_depths.'.join([str(s) for s in sel_sample])
		gt_quals = 'gt_quals.' + ', gt_quals.'.join([str(s) for s in sel_sample])
		gt_copy_numbers = 'gt_copy_numbers.' + ', gt_copy_numbers.'.join([str(s) for s in sel_sample])
		gt_phred_ll_homref = 'gt_phred_ll_homref.' + ', gt_phred_ll_homref.'.join([str(s) for s in sel_sample])
		gt_phred_ll_het = 'gt_phred_ll_het.' + ', gt_phred_ll_het.'.join([str(s) for s in sel_sample])
		gt_phred_ll_homalt = 'gt_phred_ll_homalt.' + ', gt_phred_ll_homalt.'.join([str(s) for s in sel_sample])

		gt_filter  = " != './.' or".join([s for s in gts.split(',')]) + " != './.' "
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, \
					""" + gts +', ' + gt_types +', ' + gt_phases +', ' +gt_depths +', ' + gt_ref_depths +', ' +gt_alt_depths +', ' + gt_quals +', ' + gt_copy_numbers+', ' + gt_phred_ll_homref+', ' +gt_phred_ll_het+', ' + gt_phred_ll_homalt+', ' + """\
					g.gene, g.ensembl_gene_id, g.synonym, g.transcript_min_start, g.transcript_max_end \
					from variants_cnv v, gene_view g \
					where (v.chrom ==  g.chrom \
					and g.transcript_min_start >= v.start \
					and g.transcript_max_end <= v.end) \
					or \
					(v.chrom ==  g.chrom \
					and g.transcript_min_start < v.start \
					and g.transcript_max_end > v.start) \
					or \
					(v.chrom ==  g.chrom \
					and g.transcript_min_start < v.end \
					and g.transcript_max_end > v.end) \
					order by v.chrom,v.start"""
	else :
		gt_filter = None
		sel_sample = _sample_name(database = args.db)
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end,\
					(gts).(*), (gt_types).(*),(gt_phases).(*),(gt_depths).(*),(gt_ref_depths).(*),(gt_alt_depths).(*),\
					(gt_quals).(*),(gt_copy_numbers).(*),(gt_phred_ll_homref).(*),(gt_phred_ll_het).(*),(gt_phred_ll_homalt).(*), \
					g.gene, g.ensembl_gene_id, g.synonym,g.transcript_min_start, g.transcript_max_end \
					from variants_cnv v, gene_view g \
					where (v.chrom ==  g.chrom \
					and g.transcript_min_start >= v.start \
					and g.transcript_max_end <= v.end) \
					or \
					(v.chrom ==  g.chrom \
					and g.transcript_min_start < v.start \
					and g.transcript_max_end > v.start) \
					or \
					(v.chrom ==  g.chrom \
					and g.transcript_min_start < v.end \
					and g.transcript_max_end > v.end) \
					order by v.chrom,v.start"""

	return query,gt_filter,sel_sample

def _query_custom(args):
	if args.sample:
		sel_sample = args.sample.split(',')
		gts = 'gts.' + ', gts.'.join([str(s) for s in sel_sample])
		gt_types = 'gt_types.' + ', gt_types.'.join([str(s) for s in sel_sample])
		gt_phases = 'gt_phases.' + ', gt_phases.'.join([str(s) for s in sel_sample])
		gt_depths = 'gt_depths.' + ', gt_depths.'.join([str(s) for s in sel_sample])
		gt_ref_depths = 'gt_ref_depths.' + ', gt_ref_depths.'.join([str(s) for s in sel_sample])
		gt_alt_depths = 'gt_alt_depths.' + ', gt_alt_depths.'.join([str(s) for s in sel_sample])
		gt_quals = 'gt_quals.' + ', gt_quals.'.join([str(s) for s in sel_sample])
		gt_copy_numbers = 'gt_copy_numbers.' + ', gt_copy_numbers.'.join([str(s) for s in sel_sample])
		gt_phred_ll_homref = 'gt_phred_ll_homref.' + ', gt_phred_ll_homref.'.join([str(s) for s in sel_sample])
		gt_phred_ll_het = 'gt_phred_ll_het.' + ', gt_phred_ll_het.'.join([str(s) for s in sel_sample])
		gt_phred_ll_homalt = 'gt_phred_ll_homalt.' + ', gt_phred_ll_homalt.'.join([str(s) for s in sel_sample])

		gt_filter  = " != './.' or".join([s for s in gts.split(',')]) + " != './.' "
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, \
					""" + gts +', ' + gt_types +', ' + gt_phases +', ' +gt_depths +', ' + gt_ref_depths +', ' +gt_alt_depths +', ' + gt_quals +', ' + gt_copy_numbers+', ' + gt_phred_ll_homref+', ' +gt_phred_ll_het+', ' + gt_phred_ll_homalt+', ' + """\
					g.gene, g.transcript_min_start, g.transcript_max_end
					from variants_cnv v, gene_custom_map g
					where (v.chrom ==  g.chrom
					and g.transcript_min_start >= v.start
					and g.transcript_max_end <= v.end)
					or
					(v.chrom ==  g.chrom
					and g.transcript_min_start < v.start
					and g.transcript_max_end > v.start)
					or
					(v.chrom ==  g.chrom
					and g.transcript_min_start < v.end
					and g.transcript_max_end > v.end)
					order by v.chrom,v.start"""
	else :
		gt_filter = None
		sel_sample = _sample_name(database = args.db)
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, \
					(gts).(*), (gt_types).(*),(gt_phases).(*),(gt_depths).(*),(gt_ref_depths).(*),(gt_alt_depths).(*),\
					(gt_quals).(*),(gt_copy_numbers).(*),(gt_phred_ll_homref).(*),(gt_phred_ll_het).(*),(gt_phred_ll_homalt).(*), \
		 			g.gene, g.transcript_min_start, g.transcript_max_end \
					from variants_cnv v, gene_custom_map g \
					where (v.chrom ==  g.chrom \
					and g.transcript_min_start >= v.start \
					and g.transcript_max_end <= v.end) \
					or \
					(v.chrom ==  g.chrom \
					and g.transcript_min_start < v.start \
					and g.transcript_max_end > v.start) \
					or \
					(v.chrom ==  g.chrom \
					and g.transcript_min_start < v.end \
					and g.transcript_max_end > v.end) \
					order by v.chrom,v.start"""

	return query,gt_filter,sel_sample



def overlap_gene_main(args):
	"""
	A function to view the overlap between gene and CNV
	"""
	query,gt_filter,sel_sample = _query(args)

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(query,gt_filter)

	gene[:] = []
	alt = np.array([np.zeros(len(sel_sample))])
	smp2idx = res.sample_to_idx
	result = []
	for row in res:
		gene.append(str(row['gene']))
		over_perc = perc_over(row)
		alt = _extract_data(alt,row,sel_sample,smp2idx)
		gts = row['gts']
		record = dict(variant_id = row["variant_id"],
					chrom = str(row["chrom"]),
					type = str(row["type"]),
					sub_type = str(row["sub_type"]),
					alt = str(row["alt"]),
					sv_length = row["sv_length"],
					start = row["start"],
					end = row["end"],
					gene = _genename(row),
					transcript_min_start = row["transcript_min_start"],
					transcript_max_end = row["transcript_max_end"],
					perc = over_perc,

					gts=pack_blob(gts),
					gt_types = pack_blob(row['gt_types']), gt_phases = pack_blob(row['gt_phases']) ,
                    gt_depths=pack_blob(row['gt_depths']), gt_ref_depths=pack_blob(row['gt_ref_depths']),
                    gt_alt_depths=pack_blob(row['gt_alt_depths']),
                    gt_quals=pack_blob(row['gt_quals']), gt_copy_numbers=pack_blob(row['gt_copy_numbers']),
					gt_phred_ll_homref=pack_blob(row['gt_phred_ll_homref']),
					gt_phred_ll_het=pack_blob(row['gt_phred_ll_het']),
					gt_phred_ll_homalt=pack_blob(row['gt_phred_ll_homalt']))

		if args.perc_max and args.perc_min:
			if float(args.perc_max) > over_perc and float(args.perc_min) < over_perc:
				result.append(record)
		elif args.perc_max:
			if float(args.perc_max) > over_perc:
				result.append(record)
		elif args.perc_min:
			if float(args.perc_min) < over_perc:
				result.append(record)
		else:
			result.append(record)

	alt = np.delete(alt,0,0)

	e = sql.create_engine(database.get_path(args.db), isolation_level=None)
	e.connect().connection.connection.text_factory = str
	metadata = sql.MetaData(bind=e)
	session = create_session(bind=e, autocommit=False, autoflush=False)
	gene_overlap_result = sql.Table('overlap_gene_result', metadata)

	c, metadata = database.get_session_metadata(args.db)
	c.execute("DROP TABLE if exists overlap_gene_result")
	database.create_gene_overlap_result(c,metadata, args)
	database.insert_overlap_gene_result(c, metadata,result)
	database.close_and_commit(c)

	if args.heatmap:
		heatmap(database=args.db,alt=alt,gene=gene,sel_sample=sel_sample)

def overlap_gene_browser(args):
	"""
	A function to view the overlap between gene and CNV
	"""
	result = []

	query,gt_filter,sel_sample = _query(args)

	res = GeminiQuery.GeminiQuery(args.db)
	res._set_gemini_browser(True)
	res.run(query,gt_filter)

	gene[:] = []
	alt = np.array([np.zeros(len(sel_sample))])

	for row in res:
		row['over_perc'] = perc_over(row)
		gene.append(str(row['gene']))
		temp = []
		for s in sel_sample:
			if row['gts.'+str(s)] != './.':
				if row['alt'] == 'DUP':
					temp.append(1)
				elif row['alt'] == 'DEL':
					temp.append(-1)
				else:
					temp.append(0)
			else:
				temp.append(0)
		alt = np.concatenate((alt,[temp]),axis=0)

		result.append(row)

	alt = np.delete(alt,0,0)

	return gene, alt

#######################################
## use a custom gene map to the join ##
# tabed: chrom, start, end, gene_name #
#######################################

def _get_gene_map(args):
	"""
	Define a custom gene map table
	"""
	c, metadata = database.get_session_metadata(args.db)
	# drop table if already exists
	c.execute("DROP TABLE if exists gene_custom_map")
	# create table
	database.create_gene_custom_table(c,metadata,args)
	#unique identifier for each entry
	i = 0
	contents = gene_map_c = []

	with open(args.gene_map,'r') as g:
		next(g)
		for line in g:
			col = line.strip().split("\t")
			table = gene_table.gene_custom_map(col)
			i += 1
			gene_map_c = [str(i),table.chrom,table.transcript_min_start,table.transcript_max_end,table.gene]
			contents.append(gene_map_c)
			if i % 1000 == 0:
				database.insert_gene_custom_map(c,metadata, contents)
				contents = []
		database.insert_gene_custom_map(c, metadata, contents)
		database.close_and_commit(c)

def overlap_custom_gene(args):
	_get_gene_map(args=args)

	query_custom,gt_filter,sel_sample = _query_custom(args)

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(query_custom,gt_filter)

	gene[:] = []
	alt = np.array([np.zeros(len(sel_sample))])

	smp2idx = res.sample_to_idx
	result = []
	for row in res:
		over_perc = perc_over(row)
		gene.append(str(row['gene']))
		alt = _extract_data(alt,row,sel_sample,smp2idx)
		gts = row['gts']
		record = dict(variant_id = row["variant_id"],
					chrom = str(row["chrom"]),
					type = str(row["type"]),
					sub_type = str(row["sub_type"]),
					alt = str(row["alt"]),
					sv_length = row["sv_length"],
					start = row["start"],
					end = row["end"],
					gene = row["gene"],
					transcript_min_start = row["transcript_min_start"],
					transcript_max_end = row["transcript_max_end"],
					perc = over_perc,
					gts=pack_blob(gts),
					gt_types = pack_blob(row['gt_types']), gt_phases = pack_blob(row['gt_phases']) ,
                    gt_depths=pack_blob(row['gt_depths']), gt_ref_depths=pack_blob(row['gt_ref_depths']),
                    gt_alt_depths=pack_blob(row['gt_alt_depths']),
                    gt_quals=pack_blob(row['gt_quals']), gt_copy_numbers=pack_blob(row['gt_copy_numbers']),
					gt_phred_ll_homref=pack_blob(row['gt_phred_ll_homref']),
					gt_phred_ll_het=pack_blob(row['gt_phred_ll_het']),
					gt_phred_ll_homalt=pack_blob(row['gt_phred_ll_homalt']))

		if args.perc_max and args.perc_min:
			if float(args.perc_max) > over_perc and float(args.perc_min) < over_perc:
				result.append(record)
		elif args.perc_max:
			if float(args.perc_max) > over_perc:
				result.append(record)
		elif args.perc_min:
			if float(args.perc_min) < over_perc:
				result.append(record)
		else:
			result.append(record)

	alt = np.delete(alt,0,0)

	e = sql.create_engine(database.get_path(args.db), isolation_level=None)
	e.connect().connection.connection.text_factory = str
	metadata = sql.MetaData(bind=e)
	session = create_session(bind=e, autocommit=False, autoflush=False)
	gene_overlap_result = sql.Table('overlap_gene_result', metadata)

	c, metadata = database.get_session_metadata(args.db)
	c.execute("DROP TABLE if exists overlap_gene_result")
	database.create_gene_overlap_result(c,metadata, args)
	database.insert_overlap_gene_result(c, metadata,result)
	database.close_and_commit(c)

	if args.heatmap:
		heatmap(database=args.db,alt=alt,gene = gene,sel_sample=sel_sample)

def overlap_custom_gene_browser(args):
	"""
	A function to view the overlap between gene and CNV
	"""

	result = []

	query_custom,gt_filter,sel_sample = _query(args)

	res = GeminiQuery.GeminiQuery(args.db)
	res._set_gemini_browser(True)
	res.run(query_custom,gt_filter)

	gene[:] = []
	alt = np.array([np.zeros(len(sel_sample))])

	for row in res:
		row['over_perc'] = perc_over(row)
		gene.append(str(row['gene']))
		temp = []
		for i,s in enumerate(sel_sample):
			if row['gts.'+str(s)] != './.':
				if row['alt'] == 'DUP':
					temp.append(1)
				elif row['alt'] == 'DEL':
					temp.append(-1)
				else:
					temp.append(0)
			else:
				temp.append(0)
		alt = np.concatenate((alt,[temp]),axis=0)
		result.append(row)

	alt = np.delete(alt,0,0)

	return gene, alt, result

def _sample_name(database):
	names = []
	query = "SELECT name FROM samples"
	name = GeminiQuery.GeminiQuery(database)
	name.run(query)
	for n in name:
		names.append(str(n))
	return names

def heatmap(database,alt,gene, sel_sample):
	import numpy as np
	import seaborn as sb; sb.set()

	alt_a = np.array(alt)

	# get the tick label font size
	fontsize_pt = plt.rcParams['ytick.labelsize']
	dpi = 72

	# comput the matrix height in points and inches
	matrix_height_pt = fontsize_pt * alt_a.shape[0]
	matrix_height_in = matrix_height_pt / dpi

	# compute the required figure height
	top_margin = 0.04  # in percentage of the figure height
	bottom_margin = 0.04 # in percentage of the figure height
	figure_height = matrix_height_in / (1 - top_margin - bottom_margin)

	# build the figure instance with the desired height
	fig, ax = plt.subplots(
	        figsize=(6,figure_height),
	        gridspec_kw=dict(top=1-top_margin, bottom=bottom_margin))


	ax = sb.heatmap(alt_a,linewidths=.2, ax=ax, yticklabels = False)

	plt.yticks(np.arange(len(gene)),gene, verticalalignment = 'bottom')
	ax.set_xticklabels(sel_sample)
	#ax.set_yticklabels(gene,rotation=0)
	cbar = ax.collections[0].colorbar
	cbar.set_ticks([-1, 0, 1])
	cbar.set_ticklabels(['DEL', 'none', 'DUP'])

	# save the figure
	name, ext = str(database).split('.')
	path_name = os.getcwd() + '/'
	plt.savefig(path_name + name +'_gene_heatmap.png',dpi=dpi,format = 'png')
	print path_name

def perc_over(row):
	start_var = float(row['start'])
	end_var = float(row['end'])
	len_var = float(row['sv_length'])
	start_gene = float(row['transcript_min_start'])
	end_gene = float(row['transcript_max_end'])
	len_gene = end_gene - start_gene

	if start_gene < start_var and end_gene > start_var:
		over_perc = round(float((end_gene - start_var)/len_gene)*100,2)
	elif end_gene > end_var and start_gene < end_var:
		over_perc = round(float((end_var - start_gene)/len_gene)*100,2)
	else:
		over_perc = round(float(100),2)
	return over_perc


def print_rec(args):
	args.query = "select * from overlap_gene_result"

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(args.query)
	print res.header
	for row in res:
		print row
