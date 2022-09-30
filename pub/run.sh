#!/usr/bin/env bash

#==========#
# Run Best #
#==========#

REFERENCE=reference/chm13.draft_v1.0.fasta.gz

best illumina/chm13.draft_v1.0.pcrfree.bam ${REFERENCE} illumina/illumina
best ont/chm13.draft_v1.0.ont_guppy_3.6.0.wm_2.01.pri.bam ${REFERENCE} ont/ont
best pacbio/chm13.draft_v1.0.hifi_20k.wm_2.01.pri.bam ${REFERENCE} pacbio/pacbio
