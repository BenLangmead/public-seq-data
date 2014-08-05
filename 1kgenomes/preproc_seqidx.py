#!/usr/bin/env python

# python preproc_seqidx.py sequence.index 1KG_P2_NA12878_GAII 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/' 'SAMPLE_NAME=NA12878,INSTRUMENT_PLATFORM=ILLUMINA,INSTRUMENT_MODEL=Illumina Genome Analyzer II,STUDY_NAME=1000Genomes Project Pilot 2'

import sys
import os
import string

seq_index = sys.argv[1]
group     = sys.argv[2]
prefix    = sys.argv[3]
filterstr = sys.argv[4]
filters   = []

colmap = { \
    'FASTQ_FILE'          : 1, \
    'MD5'                 : 2, \
    'RUN_ID'              : 3, \
    'STUDY_ID'            : 4, \
    'STUDY_NAME'          : 5, \
    'CENTER_NAME'         : 6, \
    'SUBMISSION_ID'       : 7, \
    'SUBMISSION_DATE'     : 8, \
    'SAMPLE_ID'           : 9, \
    'SAMPLE_NAME'         : 10, \
    'POPULATION'          : 11, \
    'EXPERIMENT_ID'       : 12, \
    'INSTRUMENT_PLATFORM' : 13, \
    'INSTRUMENT_MODEL'    : 14, \
    'LIBRARY_NAME'        : 15, \
    'RUN_NAME'            : 16, \
    'RUN_BLOCK_NAME'      : 17, \
    'INSERT_SIZE'         : 18, \
    'LIBRARY_LAYOUT'      : 19, \
    'PAIRED_FASTQ'        : 20, \
    'WITHDRAWN'           : 21, \
    'WITHDRAWN_DATE'      : 22, \
    'COMMENT'             : 23, \
    'READ_COUNT'          : 24, \
    'BASE_COUNT'          : 25, \
    'ANALYSIS_GROUP'      : 26  \
}

# Parse filters
if len(filterstr) > 0:
    print >> sys.stderr, "Parsing filter string '%s'" % filterstr
    for tok in string.split(filterstr, ','):
        kv = string.split(tok, '=')
        assert len(kv) == 2
        assert kv[0] in colmap
        filters.append((colmap[kv[0]]-1, kv[1]))
    print >> sys.stderr, "Parsed %d filters" % len(filters)

# Generate manifest from sequence.index file
fh = open(seq_index, 'r')
head = fh.readline()
while True:
    ln = fh.readline().rstrip()
    if len(ln) == 0: break
    ts = string.split(ln, '\t')
    passed = True
    for fl in filters:
        col, val = fl[0], fl[1]
        if ts[col] != val:
            passed = False
            break
    if not passed: continue
    url = prefix + ts[0]
    name = url
    name = string.split(name, '/')[-1]
    name = string.split(name, '.')[0]
    # TODO: handle paired-end
    print "\t".join([group, name, url])

fh.close()
