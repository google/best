#!/usr/bin/env bash
# Copyright (c) 2022 Google LLC
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


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
