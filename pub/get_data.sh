#!/usr/bin/env bash

#===================#
# Download Datasets #
#===================#

mkdir reference illumina ont pacbio

# Download the CHM13 v1.0 draft
wget -O reference/chm13.draft_v1.0.fasta.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.0.fasta.gz

function download_bam {
  SUBSAMPLE=${1}
  DIRECTORY=${2}
  URL=${3}
  curl ${URL} | \
  samtools view -s ${SUBSAMPLE} -bh > ${DIRECTORY}/$(basename ${URL})
}

# Illumina
download_bam 0.1 illumina https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.pcrfree.bam

# ONT
download_bam 0.1 ont https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.ont_guppy_3.6.0.wm_2.01.pri.bam

# PacBio HiFi
download_bam 0.1 pacbio https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.hifi_20k.wm_2.01.pri.bam
