use clap::Parser;

use rayon::prelude::*;

use noodles::bam;
use noodles::fasta;

use fxhash::FxHashMap;

use std::fs::File;
use std::path::PathBuf;
use std::io::{BufReader, Write};
use std::sync::Mutex;
use std::time::Instant;

mod stats;
use stats::*;
mod summary;
use summary::*;

const PER_ALN_STATS_NAME: &str = "per_aln_stats.csv";
const YIELD_STATS_NAME: &str = "summary_yield_stats.csv";
const IDENTITY_STATS_NAME: &str = "summary_identity_stats.csv";

fn run(input_path: String, reference_path: String, stats_prefix: String) {
    // read reference sequences from fasta file
    let mut ref_reader = fasta::Reader::new(BufReader::new(File::open(reference_path).unwrap()));
    let reference_seqs: FxHashMap<String, fasta::Record> = ref_reader
        .records()
        .map(|r| r.unwrap())
        .map(|r| (r.name().to_string(), r))
        .collect();

    // read bam file
    let mut reader = bam::Reader::new(File::open(input_path).unwrap());
    reader.read_header().unwrap();
    let references = reader.read_reference_sequences().unwrap();

    std::fs::create_dir_all(&stats_prefix).unwrap();

    // create per alignment stats writer that is shared between threads
    let mut aln_stats_path = PathBuf::from(&stats_prefix);
    aln_stats_path.push(PER_ALN_STATS_NAME);
    let mut aln_stats_writer = File::create(&aln_stats_path).unwrap();
    write!(aln_stats_writer, "{}\n", AlnStats::header()).unwrap();
    let aln_stats_writer = Mutex::new(aln_stats_writer);

    let summary_yield = Mutex::new(YieldSummary::new());
    let summary_identity = Mutex::new(IdentitySummary::new());

    // lazily read records to shift parsing work to individual threads
    reader
        .lazy_records()
        .par_bridge()
        .map(|r| r.unwrap())
        .for_each(|record| {
            let stats = AlnStats::from_record(&references, &reference_seqs, &record);

            if let Some(stats) = stats {
                summary_yield.lock().unwrap().update(&stats);
                summary_identity.lock().unwrap().update(&stats);

                let mut writer = aln_stats_writer.lock().unwrap();
                write!(writer, "{}\n", stats.to_csv()).unwrap();
            }
        });

    let summary_yield = summary_yield.into_inner().unwrap();
    let mut summary_yield_path = PathBuf::from(&stats_prefix);
    summary_yield_path.push(YIELD_STATS_NAME);
    let mut summary_yield_writer = File::create(&summary_yield_path).unwrap();
    write!(summary_yield_writer, "{}", summary_yield).unwrap();

    let summary_identity = summary_identity.into_inner().unwrap();
    let mut summary_identity_path = PathBuf::from(&stats_prefix);
    summary_identity_path.push(IDENTITY_STATS_NAME);
    let mut summary_identity_writer = File::create(&summary_identity_path).unwrap();
    write!(summary_identity_writer, "{}", summary_identity).unwrap();
}

fn main() {
    let args = Args::parse();
    let start_time = Instant::now();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    run(args.input, args.reference, args.stats_prefix);

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

    /// Output directory with statistics.
    stats_prefix: String,

    /// Number of threads.
    #[clap(short, long, default_value_t = 0usize)]
    threads: usize,
}
