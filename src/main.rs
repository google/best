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
mod summary;
use summary::*;
mod bed;
use bed::*;

const PER_ALN_STATS_NAME: &str = "per_aln_stats.csv";
const YIELD_STATS_NAME: &str = "summary_yield_stats.csv";
const IDENTITY_STATS_NAME: &str = "summary_identity_stats.csv";
const FEATURE_STATS_NAME: &str = "summary_feature_stats.csv";

fn run(
    input_path: String,
    reference_path: String,
    stats_prefix: String,
    intervals_path: Option<String>,
) {
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

    let intervals = intervals_path.map(|p| Intervals::new(&p));

    // create per alignment stats writer that is shared between threads
    let aln_stats_path = format!("{}.{}", stats_prefix, PER_ALN_STATS_NAME);
    let mut aln_stats_writer = File::create(&aln_stats_path).unwrap();
    write!(aln_stats_writer, "prefix,{}\n", AlnStats::header()).unwrap();
    let aln_stats_writer = Mutex::new(aln_stats_writer);

    let summary_yield = Mutex::new(YieldSummary::new(stats_prefix.clone()));
    let summary_identity = Mutex::new(IdentitySummary::new(stats_prefix.clone()));
    let summary_features = Mutex::new(FeatureSummary::new(stats_prefix.clone()));

    // lazily read records to shift parsing work to individual threads
    reader
        .lazy_records()
        .par_bridge()
        .map(|r| r.unwrap())
        .for_each(|record| {
            let aln_ref = references[record.reference_sequence_id().unwrap().unwrap()]
                .name()
                .as_str();
            // convert to one-indexed [aln_start, aln_end)
            let aln_start = usize::from(record.alignment_start().unwrap().unwrap());
            let aln_end = aln_start + record.cigar().unwrap().alignment_span();
            let overlap_intervals = intervals
                .as_ref()
                .map(|i| i.find(aln_ref, aln_start, aln_end))
                .unwrap_or_else(|| Vec::new());

            let stats =
                AlnStats::from_record(&references, &reference_seqs, &record, &overlap_intervals);

            if let Some(stats) = stats {
                summary_yield.lock().unwrap().update(&stats);
                summary_identity.lock().unwrap().update(&stats);
                if intervals.is_some() {
                    summary_features.lock().unwrap().update(&stats);
                }

                let mut writer = aln_stats_writer.lock().unwrap();
                write!(writer, "{},{}\n", stats_prefix, stats.to_csv()).unwrap();
            }
        });

    let summary_yield = summary_yield.into_inner().unwrap();
    let summary_yield_path = format!("{}.{}", stats_prefix, YIELD_STATS_NAME);
    let mut summary_yield_writer = File::create(&summary_yield_path).unwrap();
    write!(summary_yield_writer, "{}", summary_yield).unwrap();

    let summary_identity = summary_identity.into_inner().unwrap();
    let summary_identity_path = format!("{}.{}", stats_prefix, IDENTITY_STATS_NAME);
    let mut summary_identity_writer = File::create(&summary_identity_path).unwrap();
    write!(summary_identity_writer, "{}", summary_identity).unwrap();

    if intervals.is_some() {
        let summary_features = summary_features.into_inner().unwrap();
        let summary_features_path = format!("{}.{}", stats_prefix, FEATURE_STATS_NAME);
        let mut summary_features_writer = File::create(&summary_features_path).unwrap();
        write!(summary_features_writer, "{}", summary_features).unwrap();
    }
}

fn main() {
    let args = Args::parse();
    let start_time = Instant::now();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    run(
        args.input,
        args.reference,
        args.stats_prefix,
        args.intervals,
    );

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

    /// Input intervals BED file.
    intervals: Option<String>,

    /// Output directory with statistics.
    stats_prefix: String,

    /// Number of threads. Will be automatically determined if this is set to 0.
    #[clap(short, long, default_value_t = 0usize)]
    threads: usize,
}
