#!/bin/sh

if [ ! -f "sequence.index" ] ; then
	INDEXES="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/alignment.index "
	INDEXES="$INDEXES ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/sequence.index"
	for i in $INDEXES ; do wget $i ; done
fi

python preproc_seqidx.py \
	sequence.index \
	1KG_P2_NA12878_GAII \
	'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/' \
	'SAMPLE_NAME=NA12878,INSTRUMENT_PLATFORM=ILLUMINA,INSTRUMENT_MODEL=Illumina Genome Analyzer II,STUDY_NAME=1000Genomes Project Pilot 2' \
	> 1KG_P2_NA12878_GAII.manifest

python preproc_seqidx.py \
	sequence.index \
	1KG_ASW_HISEQ \
	'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/' \
	'INSTRUMENT_PLATFORM=ILLUMINA,INSTRUMENT_MODEL=Illumina HiSeq 2000,STUDY_NAME=1000 Genomes CEPH (Utah residents with ancestry from Northern and Western Europe) population sequencing' \
	> 1KG_ASW_HISEQ.manifest
