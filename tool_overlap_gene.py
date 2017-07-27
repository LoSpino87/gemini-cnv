#A tool for valutate overlapping of genome
#!/usr/bin/env python

# Python imports
import os.path
import os
import re
import sys
import matplotlib
import matplotlib.pyplot as plt

# Gemini imports
import GeminiQuery

gene = []
alt = []
samples = []
result = []

def run(parser, args):
	if os.path.exists(args.db):
		overlap_gene_main(args)
		samples = sample_name(database = args.db)

def overlap_gene_main(args):
	import numpy as np
	"""
	A function to view the overlap between gene and CNV
	"""
	args.query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.*
				from variants_cnv v, gene_view g
				where g.chrom == v.chrom
				and g.transcript_min_start >= v.start
				and g.transcript_max_end <= v.end
				order by g.gene"""

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(args.query)

	for row in res:
		gene.append(str(row['gene']))
		if row['alt'] == 'DUP':
			alt.append(1)
		else: alt.append(-1)
		print row

	if args.heatmap:
		heatmap(database=args.db)
		name, ext = str(args.db).split('.')
		path_name = os.getcwd() + '/' + name + '/'
		picture = path_name + name + '_overlap_gene.png'


def overlap_gene_browser(database):

	"""
	A function to view the overlap between gene and CNV
	"""
	query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.*
				from variants_cnv v, gene_view g
				where g.chrom == v.chrom
				and g.transcript_min_start >= v.start
				and g.transcript_max_end <= v.end
				order by g.gene"""
	res = GeminiQuery.GeminiQuery(database)
	res._set_gemini_browser(True)
	res.run(query)

	for row in res:
		gene.append(str(row['gene']))
		if row['alt'] == 'DUP':
			alt.append(1)
		else: alt.append(-1)
		result.append(row)

	return result

def sample_name(database):
	names = []
	query = "SELECT name FROM samples"
	name = GeminiQuery.GeminiQuery(database)
	name.run(query)
	for n in name:
		names.append(n)
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
	path_name = os.getcwd() + '/' + name + '/'
	plt.savefig(path_name + name +'_overlap_gene.png')
