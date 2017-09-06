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

# Gemini imports
import GeminiQuery
import gene_table
import database

gene = list()
alt = list()
samples = []

def run(parser, args):
	if os.path.exists(args.db):
		samples = sample_name(database = args.db)

		if args.gene_map:
			overlap_custom_gene(args)
		else:
			overlap_gene_main(args)

def overlap_gene_main(args):
	"""
	A function to view the overlap between gene and CNV
	"""
	gene[:] = []
	alt[:] = []
	if args.sample:
		gt_col = 'gts.' + ', gts.'.join([s for s in str(args.sample).split(',')])
		gt_filter  = " != './.' or".join([s for s in gt_col.split(',')]) + " != './.' "
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, """ + gt_col + """, g.gene, g.ensembl_gene_id, g.synonym
					from variants_cnv v, gene_view g
					where g.chrom == v.chrom
					and g.transcript_min_start >= v.start
					and g.transcript_max_end <= v.end
					order by g.gene"""
	else :
		gt_filter = None
		query = query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.gene, g.ensembl_gene_id, g.synonym
					from variants_cnv v, gene_view g
					where g.chrom == v.chrom
					and g.transcript_min_start >= v.start
					and g.transcript_max_end <= v.end
					order by g.gene"""

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(query,gt_filter)

	for row in res:
		gene.append(str(row['gene']))
		if row['alt'] == 'DUP':
			alt.append(1)
		elif row['alt'] == 'DEL':
			alt.append(-1)
		else:
			alt.append(0)
		print row

	if args.heatmap:
		heatmap(database=args.db)


def overlap_gene_browser(args):
	"""
	A function to view the overlap between gene and CNV
	"""
	result = []
	gene[:] = []
	alt[:] = []

	if args.sample:
		gt_col = 'gts.' + ', gts.'.join([s for s in str(args.sample).split(',')])
		gt_filter  = " != './.' or".join([s for s in gt_col.split(',')]) + " != './.' "
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, """ + gt_col + """, g.gene, g.ensembl_gene_id, g.synonym
					from variants_cnv v, gene_view g
					where g.chrom == v.chrom
					and g.transcript_min_start >= v.start
					and g.transcript_max_end <= v.end
					order by g.gene"""
	else :
		gt_filter = None
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.gene, g.ensembl_gene_id, g.synonym
					from variants_cnv v, gene_view g
					where g.chrom == v.chrom
					and g.transcript_min_start >= v.start
					and g.transcript_max_end <= v.end
					order by g.gene"""

	res = GeminiQuery.GeminiQuery(args.db)
	res._set_gemini_browser(True)
	res.run(query,gt_filter)

	for row in res:
		gene.append(str(row['gene']))
		if row['alt'] == 'DUP':
			alt.append(1)
		elif row['alt'] == 'DEL':
			alt.append(-1)
		else:
			alt.append(0)
		result.append(row)

	return result

#######################################
## use a custom gene map to the join ##
# tabed: chrom, start, end, gene_name #
#######################################

def overlap_custom_gene(args):
	get_gene_map(args=args)

	gene[:] = []
	alt[:] = []

	if args.sample:
		gt_col = 'gts.' + ', gts.'.join([s for s in str(args.sample).split(',')])
		gt_filter  = " != './.' or".join([s for s in gt_col.split(',')]) + " != './.' "
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, """ + gt_col + """, g.gene, g.ensembl_gene_id, g.synonym
					from variants_cnv v, gene_custom_map g
					where g.chrom == v.chrom
					and g.start >= v.start
					and g.end <= v.end
					order by g.gene_name"""
	else :
		gt_filter = None
		query_custom = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.gene_name
					from variants_cnv v, gene_custom_map g
					where g.chrom == v.chrom
					and g.start >= v.start
					and g.end <= v.end
					order by g.gene_name"""

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(query_custom,gt_filter)

	for row in res:
		gene.append(str(row['gene_name']))
		if row['alt'] == 'DUP':
			alt.append(1)
		elif row['alt'] == 'DEL':
			alt.append(-1)
		else:
			alt.append(0)
		print row

	if args.heatmap:
		heatmap(database=args.db)


def overlap_custom_gene_browser(args):
	"""
	A function to view the overlap between gene and CNV
	"""

	gene[:] = []
	alt[:] = []
	result = []

	if args.sample:
		gt_col = 'gts.' + ', gts.'.join([s for s in str(args.sample).split(',')])
		gt_filter  = " != './.' or".join([s for s in gt_col.split(',')]) + " != './.' "
		query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, """ + gt_col + """, g.gene, g.ensembl_gene_id, g.synonym
					from variants_cnv v, gene_custom_map g
					where g.chrom == v.chrom
					and g.start >= v.start
					and g.end <= v.end
					order by g.gene_name"""
	else :
		gt_filter = None
		query_custom = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.gene_name
					from variants_cnv v, gene_custom_map g
					where g.chrom == v.chrom
					and g.start >= v.start
					and g.end <= v.end
					order by g.gene_name"""

	res = GeminiQuery.GeminiQuery(args.db)
	res._set_gemini_browser(True)
	res.run(query_custom,gt_filter)

	for row in res:
		gene.append(str(row['gene_name']))
		if row['alt'] == 'DUP':
			alt.append(1)
		elif row['alt'] == 'DEL':
			alt.append(-1)
		else:
			alt.append(0)
		result.append(row)

	return result

def get_gene_map(args):
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
			gene_map_c = [str(i),table.chrom,table.start,table.end,table.gene_name]
			contents.append(gene_map_c)
			if i % 1000 == 0:
				database.insert_gene_custom_map(c,metadata, contents)
				contents = []
		database.insert_gene_custom_map(c, metadata, contents)
		database.close_and_commit(c)


def sample_name(database):
	names = []
	query = "SELECT name FROM samples"
	name = GeminiQuery.GeminiQuery(database)
	name.run(query)
	for n in name:
		names.append(str(n))
	return names

def heatmap(database):
	import numpy as np
	import seaborn as sb; sb.set()

	alt_a = np.array(alt)
	alt_a_T = alt_a[np.newaxis, :].T

	# get the tick label font size
	fontsize_pt = plt.rcParams['ytick.labelsize']
	dpi = 72.27

	# comput the matrix height in points and inches
	matrix_height_pt = fontsize_pt * alt_a_T.shape[0]
	matrix_height_in = matrix_height_pt / dpi

	# compute the required figure height
	top_margin = 0.04  # in percentage of the figure height
	bottom_margin = 0.04 # in percentage of the figure height
	figure_height = matrix_height_in / (1 - top_margin - bottom_margin)


	# build the figure instance with the desired height
	fig, ax = plt.subplots(
	        figsize=(6,figure_height),
	        gridspec_kw=dict(top=1-top_margin, bottom=bottom_margin))


	ax = sb.heatmap(alt_a_T,linewidths=.2, ax=ax)
	ax.set_xticklabels(samples)
	ax.set_yticklabels(gene,rotation=0)
	cbar = ax.collections[0].colorbar
	cbar.set_ticks([-1, 0, 1])
	cbar.set_ticklabels(['DEL', 'none', 'DUP'])

	# save the figure
	name, ext = str(database).split('.')
	path_name = os.getcwd() + '/'
	print path_name
	plt.savefig(path_name + name +'_overlap_gene.png')
