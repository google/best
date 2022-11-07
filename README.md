# best
Bam Error Stats Tool (best): analysis of error types in aligned reads.

`best` is used to assess the quality of reads after aligning them to a
reference assembly.

## Features
* Collect overall and per alignment stats
* Distribution of indel lengths
* Yield at different empirical Q-value thresholds
* Bin per read stats to easily examine the distribution of errors for certain
  types of reads
* Stats for regions specified by intervals (BED file, homopolymer regions,
  windows etc.)
* Multithreading for speed

## Usage
The [`best` Usage Guide](Usage.md) gives an overview of how to use `best`.

## Development
### Running
1. Install [Rust](https://www.rust-lang.org/tools/install).
2. Run `cargo build --release`.
3. Run `cargo run --release -- input.bam reference.fasta prefix/path` or
`target/release/best input.bam reference.fasta prefix/path`.

This will generate stats files with the `prefix/path` prefix.

The built binary is located at `target/release/best`.

### Formatting
```
cargo fmt
```

### Comparing
Remember to pass the `-t 1` option to ensure that only one thread is used for
testing. Best generally tries to ensure the order of outputs is deterministic
with multiple threads, but the order of per-alignment stats is arbitrary unless
only one thread is used.
