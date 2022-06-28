use clap::Parser;

use rayon::prelude::*;

use noodles::bam;
use noodles::fasta;

use fxhash::FxHashMap;

use std::fs::File;
use std::io::{BufReader, Write};
use std::sync::Mutex;
use std::time::Instant;

mod stats;
use stats::*;

fn run(input_path: String, reference_path: String, aln_stats_path: String) {
    let mut ref_reader = fasta::Reader::new(BufReader::new(File::open(reference_path).unwrap()));
    let reference_seqs: FxHashMap<String, fasta::Record> = ref_reader
        .records()
        .map(|r| r.unwrap())
        .map(|r| (r.name().to_string(), r))
        .collect();

    let mut reader = bam::Reader::new(File::open(input_path).unwrap());
    reader.read_header().unwrap();
    let references = reader.read_reference_sequences().unwrap();

    let mut aln_stats_writer = File::create(aln_stats_path).unwrap();
    write!(aln_stats_writer, "{}\n", AlnStats::header()).unwrap();
    let aln_stats_writer = Mutex::new(aln_stats_writer);

    reader
        .lazy_records()
        .par_bridge()
        .map(|r| r.unwrap())
        .for_each(|record| {
            let stats = AlnStats::from_record(&references, &reference_seqs, &record);

            if let Some(stats) = stats {
                let mut writer = aln_stats_writer.lock().unwrap();
                write!(writer, "{}\n", stats.to_csv()).unwrap();
            }
        });
}

fn main() {
    let args = Args::parse();
    let start_time = Instant::now();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    run(args.input, args.reference, args.aln_stats);

    let duration = start_time.elapsed();
    println!("Run time (s): {}", duration.as_secs());
}

#[derive(Parser)]
#[clap(author, version, about)]
struct Args {
    /// Input BAM file.
    input: String,

    /// Input reference FASTA file.
    reference: String,

    /// Output CSV file with per alignment statistics.
    aln_stats: String,

    /// Number of threads.
    #[clap(short, long, default_value_t = 0usize)]
    threads: usize,
}
