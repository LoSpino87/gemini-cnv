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
	"""
	A function to view the overlap between gene and CNV
	"""
	args.query = """select v.variant_id, v.chrom, v.type, v.sub_type, v.alt, v.sv_length, v.start, v.end, g.*
                    from variants_cnv v, gene_summary g
                    where g.chrom == v.chrom
                    and g.transcript_min_start >= v.start
                    and g.transcript_max_end <= v.end"""

	res = GeminiQuery.GeminiQuery(args.db)
	res.run(args.query)

	for row in res:
		print row
