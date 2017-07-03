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
				from variants_cnv v, gene_summary g
				where g.chrom == v.chrom
				and g.transcript_min_start >= v.start
				and g.transcript_max_end <= v.end
				and g.ensembl_gene_id IN (SELECT d.ensembl_gene_id from gene_summary d group by d.ensembl_gene_id having count(*)=1
				UNION
				SELECT d.ensembl_gene_id from gene_summary d where d.is_hgnc=='1' group by d.ensembl_gene_id having count(*)>1 )
				group by g.ensembl_gene_id"""

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(args.query)

	sample = 'NA07000'
	gene = []
	alt = []

	for row in res:
		gene.append(row['gene'])
		if row['alt'] == 'DUP':
			alt.append(1)
		else: alt.append(0)
		print row
	#heatmap(gene,alt,sample)

def overlap_gene_browser(database):

	"""
	A function to view the overlap between gene and CNV
	"""
	query = """SELECT v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.*
				from variants_cnv v, gene_summary g
				where g.chrom == v.chrom
				and g.transcript_min_start >= v.start
				and g.transcript_max_end <= v.end
				and g.ensembl_gene_id IN (SELECT d.ensembl_gene_id from gene_summary d group by d.ensembl_gene_id having count(*)=1
				UNION
				SELECT d.ensembl_gene_id from gene_summary d where d.is_hgnc=='1' group by d.ensembl_gene_id having count(*)>1 )
				group by g.ensembl_gene_id"""
	res = GeminiQuery.GeminiQuery(database)
	res._set_gemini_browser(True)
	res.run(query)

	return res


def heatmap(gene,alt,sample):
	import matplotlib.pyplot as plt
	import matplotlib
	import numpy as np

	alt_a = np.array(alt)
	alt_a_T = alt_a[np.newaxis, :].T

	# Plot it out
	fig, ax = plt.subplots()
	heatmap = ax.pcolor(alt_a_T, cmap=plt.cm.Blues)

	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(alt_a_T[0]) + 0.5, minor=False)

	ax.set_frame_on(False)
	ax.xaxis.tick_top()

	# note I could have used nba_sort.columns but made "labels" instead
	ax.set_xticklabels(sample, minor=False)
	ax.set_yticklabels(gene, minor=False)

	# rotate the
	#plt.yticks(rotation=90)

	ax.grid(False)

	# Turn off all the ticks
	ax = plt.gca()

	for t in ax.xaxis.get_major_ticks():
	    t.tick1On = False
	    t.tick2On = False
	for t in ax.yaxis.get_major_ticks():
	    t.tick1On = False
	    t.tick2On = False

	plt.show(heatmap)
