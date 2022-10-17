# `best` Usage Guide

This guide will give a general overview of how to use `best`.

## Example Analysis
Let's say we have aligned reads in `aln.bam` and the reference assembly
`ref.fasta.gz`. We want to collect statistics on the types of
errors that occur in the alignments. Then, we can run `best` like
```
best -t 4 aln.bam ref.fasta.gz output
```
This will use 4 threads to collect stats and generate the following files:
```
output.per_aln_stats.csv.gz
output.summary_cigar_stats.csv
output.summary_identity_stats.csv
output.summary_yield_stats.csv
```
The per-alignment stats file is gzipped to save space. The
`output.summary_cigar_stats.csv` file contains the distribution of lengths of
consecutive insertions and deletions. The `output.summary_identity_stats.csv`
file contains some general stats on the error rates across all alignments. The
`output.summary_yield_stats.csv` contains the yield (number of reads/bases)
above certain quality thresholds.

We can collect even more data in one run of `best`. If we run
```
best -t 4 --intervals-hp --bin-types q_len:1000 gc_contant:0.05 -- aln.bam ref.fasta.gz output
```
then we will get two more files:
```
output.summary_feature_stats.csv
output.summary_bin_stats.csv
```
The `output.summary_feature_stats.csv` file contains stats stratified by
intervals. In this case, we use `--intervals-hp` to indicate that the intervals
are homopolymer regions, so this will produce the error types at homopolymer
regions of different lengths. The `output.summary_bin_stats.csv` file will
contain the error types for reads binned by both the read (query) length (`q_len`,
bin by increments of 1000bp) and the GC content (`gc_content`, bin by increments
of 0.05).

It is also possible to use bed files as custom intervals. These files can have
three columns (all intervals will have the same feature) or four columns, where
the last column indicates the feature. The feature stats are aggregated across
all bed intervals with the same feature.

## Help Message:
```
best 0.1.0
Daniel Liu
Bam Error Stats Tool (best): analysis of error types in aligned reads.

USAGE:
    best [OPTIONS] <INPUT> <REFERENCE> <STATS_PREFIX>

ARGS:
    <INPUT>
            Input BAM file

    <REFERENCE>
            Input reference FASTA file. Can be gzipped

    <STATS_PREFIX>
            Prefix for output files that contain statistics

OPTIONS:
    -b, --bin-types <BIN_TYPES>...
            Types of bins to use for per alignment stats.

            Each bin should be of the format <bin_type>:<step_size>.

            Supported bin types: q_len (read sequence length), subread_passes, mapq, mean_qual,
            gc_content, concordance_qv (phred scale Q-value)

    -h, --help
            Print help information

        --intervals-bed <INTERVALS_BED>...
            Use intervals from a BED file.

            The BED file should have the columns chrom, start, stop, and feature. The feature column
            is optional. It allows stats to be gathered separately for different types of intervals.

        --intervals-border <INTERVALS_BORDER>...
            Use fixed-width window border regions as intervals

        --intervals-hp
            Use homopolymer regions as intervals

        --intervals-match <INTERVALS_MATCH>...
            Use regions that match any of the specified subsequences as intervals

        --intervals-window <INTERVALS_WINDOW>...
            Use fixed-width windows as intervals

        --intervals-window-pos <INTERVALS_WINDOW_POS>...
            Use fixed-width windows with positions as intervals

    -n, --name-column <NAME_COLUMN>
            Add column with a specific name in CSV outputs

        --no-per-aln-stats
            Turn off outputting per alignment stats

    -t, --threads <THREADS>
            Number of threads. Will be automatically determined if this is set to 0

            [default: 0]

    -V, --version
            Print version information
```
