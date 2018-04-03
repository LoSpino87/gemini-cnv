#!/usr/bin/env python

import sqlalchemy as sql
from sqlalchemy.orm import mapper, create_session
import sqlalchemy

import GeminiQuery
from gemini.gim import (AutoDom, AutoRec, DeNovo, MendelViolations, CompoundHet)

database = None

def export_vcf(from_table, database, out_file):
    # write vcf header into a file
    tmp = open(out_file, 'w')

    handle = GeminiQuery.GeminiQuery(database)
    handle._set_gemini_browser(True)
    handle.run('select distinct * from vcf_header')

    for h in handle:
        lines = str(h).split('\n')
        for l in range(0,len(lines)-1):
            tmp.write(lines[l]+'\n')
        tmp.write('##Filtered with GEMINI_CNV.\n')
        tmp.write(lines[-1]+'\n')



    tmp.close
