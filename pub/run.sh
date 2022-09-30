#!/usr/bin/env bash

#===================#
# Download Datasets #
#===================#

mkdir reference illumina ont pacbio

# Download the CHM13 v1.0 draft
wget -P reference wget -P reference https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.0.fasta.gz

REFERENCE=reference/chm13.draft_v1.0.fasta.gz

function download_bam {
  SUBSAMPLE=${1}
  DIRECTORY=${2}
  URL=${3}
  curl -P ${DIRECTORY} ${URL} | \
  samtools view -s ${SUBSAMPLE} -bh > ${DIRECTORY}/$(basename ${URL})
}

# Illumina
download_bam 0.1 illumina https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.pcrfree.bam

# ONT
download_bam 0.1 ont https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.ont_guppy_3.6.0.wm_2.01.pri.bam

# PacBio HiFi
download_bam 0.1 pacbio https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.hifi_20k.wm_2.01.pri.bam

#==========#
# Run Best #
#==========#

best illumina/chm13.draft_v1.0.pcrfree.bam ${REFERENCE} illumina/illumina
best ont/chm13.draft_v1.0.ont_guppy_3.6.0.wm_2.01.pri.bam ${REFERENCE} ont/ont
best pacbio/chm13.draft_v1.0.hifi_20k.wm_2.01.pri.bam ${REFERENCE} pacbio/pacbio
