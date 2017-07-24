#A tool for valutate overlapping of genome
#!/usr/bin/env python

# Python imports
import os.path
import re
import sys

# Gemini imports
import GeminiQuery

def run(parser, args):
	if os.path.exists(args.db):
		overlap_gene_main(args)

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
				order by g.ensembl_gene_id"""

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(args.query)

	samples = sample_name(args)
	gene = []
	alt = []

	for row in res:
		gene.append(str(row['gene']))
		if row['alt'] == 'DUP':
			alt.append(1)
		else: alt.append(-1)
		print row
	heatmap(gene,alt,samples)


def overlap_gene_browser(database):

	"""
	A function to view the overlap between gene and CNV
	"""
	query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.*
				from variants_cnv v, gene_view g
				where g.chrom == v.chrom
				and g.transcript_min_start >= v.start
				and g.transcript_max_end <= v.end
				order by g.ensembl_gene_id"""
	res = GeminiQuery.GeminiQuery(database)
	res._set_gemini_browser(True)
	res.run(query)

	return res

def sample_name(args):
	names = []
	args.query = "SELECT name FROM samples"
	name = GeminiQuery.GeminiQuery(args.db)
	name.run(args.query)
	for n in name:
		names.append(n)
	return names

def heatmap(gene,alt,sample):
	import matplotlib.pyplot as plt
	import matplotlib
	import numpy as np
	import seaborn as sb; sb.set()

	alt_a = np.array(alt)
	alt_a_T = alt_a[np.newaxis, :].T

	ax = sb.heatmap(alt_a_T,linewidths=.2)
	ax.set_xticklabels(sample)
	ax.set_yticklabels(gene,rotation=0)
	plt.show(ax)
