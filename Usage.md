# `best` Usage Guide

This guide will give a general overview of how to use `best`.

## Help Message:
```
best 0.1.0
Daniel Liu
Bam Error STats (BEST): analysis of error types in aligned reads.

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
            Use intervals from a BED file

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
