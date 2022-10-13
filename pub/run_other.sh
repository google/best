#!/usr/bin/env bash

REFERENCE=reference/chm13.draft_v1.0.fasta
BAM_FILE=pacbio/chm13.draft_v1.0.hifi_20k.wm_2.01.pri.bam

BAM_CONCORDANCE=~/hg002-ccs/concordance/bamConcordance
echo "bamConcordance"
time ${BAM_CONCORDANCE} ${REFERENCE} ${BAM_FILE} pacbio/bamConcordance.csv

echo "pomoxis"
time stats_from_bam ${BAM_FILE} -o pacbio/pomoxis.tsv -s pacbio/pomoxis_summary.txt -t 4
