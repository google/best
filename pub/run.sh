#!/usr/bin/env bash

#==========#
# Run Best #
#==========#

REFERENCE=reference/chm13.draft_v1.0.fasta.gz
ARGS="-t 4 --intervals-hp --bin-types gc_content:0.05"

best ${ARGS} --bin-types q_len:10 -- illumina/chm13.draft_v1.0.pcrfree.bam ${REFERENCE} illumina/illumina
best ${ARGS} --bin-types q_len:10000 -- ont/chm13.draft_v1.0.ont_guppy_3.6.0.wm_2.01.pri.bam ${REFERENCE} ont/ont
best ${ARGS} --bin-types q_len:1000 -- pacbio/chm13.draft_v1.0.hifi_20k.wm_2.01.pri.bam ${REFERENCE} pacbio/pacbio
