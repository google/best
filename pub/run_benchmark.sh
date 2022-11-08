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
