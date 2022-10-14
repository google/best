#!/usr/bin/env bash

REFERENCE=reference/chm13.draft_v1.0.fasta
BAM_FILE=pacbio/chm13.draft_v1.0.hifi_20k.wm_2.01.pri.bam

BAM_CONCORDANCE=~/hg002-ccs/concordance/bamConcordance
echo "bamConcordance"
time ${BAM_CONCORDANCE} ${REFERENCE} ${BAM_FILE} pacbio/bamConcordance.csv

#echo "pomoxis"
#time assess_homopolymers count ${BAM_FILE} -o pacbio/pomoxis -t 4

echo "best 1 thread"
time best -t 1 ${BAM_FILE} ${REFERENCE} pacbio/best_timing

echo "best 4 thread"
time best -t 4 ${BAM_FILE} ${REFERENCE} pacbio/best_timing
