# BAM Error Stats
Analysis of error types in BAM files.

## Development
### Running
1. Install [Rust](https://www.rust-lang.org/tools/install).
2. Run `cargo build --release`.
3. Run `cargo run --release -- input.bam reference.fasta bam-error-stats` or
`target/release/bam-error-stats input.bam reference.fasta bam-error-stats`.

This will generate stats files in the `bam-error-stats` directory.

The built binary is located at `target/release/bam-error-stats`.

### Formatting
```
cargo fmt
```

### Comparing
Remember to pass the `-t 1` option to ensure that only one thread is used for
testing.
