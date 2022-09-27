#!/usr/bin/env bash

#===================#
# Download Datasets #
#===================#

mkdir reference illumina ont pacbio

# Download the CHM13 v1.0 draft
wget -P reference wget -P reference https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.0.fasta.gz
gzip -dc reference/chm13.draft_v1.0.fasta.gz > reference/chm13.draft_v1.0.fasta

N_LINES=100000
REFERENCE=reference/chm13.draft_v1.0.fasta

function download_bam {
  DIRECTORY=${1}
  URL=${2}
  curl -P ${DIRECTORY} ${URL} | \
  samtools view -h | \
  head -n ${N_LINES} | \
  samtools view -bh > ${DIRECTORY}/$(basename ${URL})
}

# Illumina
download_bam illumina https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.pcrfree.bam

# ONT
download_bam ont https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.ont_guppy_3.6.0.wm_2.01.pri.bam

# PacBio HiFi
download_bam pacbio https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/alignments/chm13.draft_v1.0.hifi_20k.wm_2.01.pri.bam

best illumina/chm13.draft_v1.0.pcrfree.bam reference./
