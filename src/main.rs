use clap::Parser;

use rayon::prelude::*;

use noodles::bam;
use noodles::fasta;

use fxhash::FxHashMap;

use std::fs::File;
use std::io::{BufReader, Write};
use std::str::FromStr;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

mod stats;
use stats::*;
mod summary;
use summary::*;
mod bed;
use bed::*;
mod intervals;
use intervals::*;

const PER_ALN_STATS_NAME: &str = "per_aln_stats.csv";
const YIELD_STATS_NAME: &str = "summary_yield_stats.csv";
const IDENTITY_STATS_NAME: &str = "summary_identity_stats.csv";
const FEATURE_STATS_NAME: &str = "summary_feature_stats.csv";
const CIGAR_STATS_NAME: &str = "summary_cigar_stats.csv";
const BIN_STATS_NAME: &str = "summary_bin_stats.csv";

fn run(
    input_path: String,
    reference_path: String,
    stats_prefix: String,
    bin_types: Option<Vec<BinType>>,
    intervals_types: Vec<IntervalsType>,
    name_column: Option<String>,
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

    // create per alignment stats writer that is shared between threads
    let aln_stats_path = format!("{}.{}", stats_prefix, PER_ALN_STATS_NAME);
    let mut aln_stats_writer = File::create(&aln_stats_path).unwrap();
    write!(
        aln_stats_writer,
        "{}{}\n",
        if name_column.is_some() { "name," } else { "" },
        AlnStats::header()
    )
    .unwrap();
    let aln_stats_writer = Mutex::new(aln_stats_writer);

    let summary_yield = Mutex::new(YieldSummary::new(name_column.clone()));
    let summary_identity = Mutex::new(IdentitySummary::new(name_column.clone()));
    let summary_features = if intervals_types.is_empty() {
        None
    } else {
        Some(Mutex::new(FeatureSummary::new(name_column.clone())))
    };
    let summary_cigars = Mutex::new(CigarLenSummary::new(name_column.clone()));
    let summary_bins = bin_types.map(|b| Mutex::new(BinSummary::new(name_column.clone(), b)));
    let total_alns = AtomicUsize::new(0);

    // lazily read records to shift parsing work to individual threads
    reader
        .lazy_records()
        .par_bridge()
        .map(|r| r.unwrap())
        .for_each(|record| {
            total_alns.fetch_add(1, Ordering::Relaxed);

            let flags = record.flags().unwrap();
            if flags.is_unmapped() || flags.is_secondary() {
                // skip
                return;
            }

            let aln_ref = references[record.reference_sequence_id().unwrap().unwrap()]
                .name()
                .as_str();
            // convert to one-indexed [aln_start, aln_end)
            let aln_start = usize::from(record.alignment_start().unwrap().unwrap());
            let aln_end = aln_start + record.cigar().unwrap().alignment_span();
            let mut intervals_vec = Vec::new();
            let mut overlap_intervals = Vec::new();
            intervals_types
                .iter()
                .for_each(|intervals_type| match intervals_type {
                    IntervalsType::Homopolymer => intervals_vec.extend(find_homopolymers(
                        reference_seqs[aln_ref].sequence(),
                        aln_start,
                        aln_end,
                    )),
                    IntervalsType::Window(win_len) => {
                        intervals_vec.extend(get_windows(aln_start, aln_end, *win_len))
                    }
                    IntervalsType::Border(win_len) => {
                        intervals_vec.extend(get_borders(aln_start, aln_end, *win_len))
                    }
                    IntervalsType::Match(seq) => intervals_vec.extend(get_matches(
                        reference_seqs[aln_ref].sequence(),
                        aln_start,
                        aln_end,
                        seq,
                    )),
                    IntervalsType::Bed(intervals) => {
                        overlap_intervals.extend(intervals.find(aln_ref, aln_start, aln_end))
                    }
                });
            overlap_intervals.extend(&intervals_vec);
            overlap_intervals.sort();

            let stats =
                AlnStats::from_record(&references, &reference_seqs, &record, &overlap_intervals);

            if let Some(stats) = stats {
                summary_yield.lock().unwrap().update(&stats);
                summary_identity.lock().unwrap().update(&stats);
                summary_features
                    .as_ref()
                    .map(|f| f.lock().unwrap().update(&stats));
                summary_cigars.lock().unwrap().update(&stats);
                summary_bins
                    .as_ref()
                    .map(|b| b.lock().unwrap().update(&stats));

                let mut writer = aln_stats_writer.lock().unwrap();
                if let Some(ref name) = name_column {
                    write!(writer, "{},{}\n", name, stats.to_csv()).unwrap();
                } else {
                    write!(writer, "{}\n", stats.to_csv()).unwrap();
                }
            }
        });

    let summary_yield = summary_yield.into_inner().unwrap();
    let summary_yield_path = format!("{}.{}", stats_prefix, YIELD_STATS_NAME);
    let mut summary_yield_writer = File::create(&summary_yield_path).unwrap();
    write!(summary_yield_writer, "{}", summary_yield).unwrap();

    let mut summary_identity = summary_identity.into_inner().unwrap();
    summary_identity.total_alns = total_alns.into_inner();
    let summary_identity_path = format!("{}.{}", stats_prefix, IDENTITY_STATS_NAME);
    let mut summary_identity_writer = File::create(&summary_identity_path).unwrap();
    write!(summary_identity_writer, "{}", summary_identity).unwrap();

    if let Some(f) = summary_features {
        let summary_features = f.into_inner().unwrap();
        let summary_features_path = format!("{}.{}", stats_prefix, FEATURE_STATS_NAME);
        let mut summary_features_writer = File::create(&summary_features_path).unwrap();
        write!(summary_features_writer, "{}", summary_features).unwrap();
    }

    let summary_cigars = summary_cigars.into_inner().unwrap();
    let summary_cigars_path = format!("{}.{}", stats_prefix, CIGAR_STATS_NAME);
    let mut summary_cigars_writer = File::create(&summary_cigars_path).unwrap();
    write!(summary_cigars_writer, "{}", summary_cigars).unwrap();

    if let Some(b) = summary_bins {
        let summary_bins = b.into_inner().unwrap();
        let summary_bins_path = format!("{}.{}", stats_prefix, BIN_STATS_NAME);
        let mut summary_bins_writer = File::create(&summary_bins_path).unwrap();
        write!(summary_bins_writer, "{}", summary_bins).unwrap();
    }
}

fn main() {
    let start_time = Instant::now();
    let args = Args::parse();

    let bin_types = args
        .bin_types
        .map(|b| b.iter().map(|s| BinType::from_str(s).unwrap()).collect());

    let mut intervals_types = Vec::new();
    if args.intervals_hp {
        intervals_types.push(IntervalsType::Homopolymer);
    }
    if let Some(paths) = args.intervals_bed {
        intervals_types.extend(paths.iter().map(|p| IntervalsType::Bed(Intervals::new(p))));
    }
    if let Some(win_lens) = args.intervals_window {
        intervals_types.extend(win_lens.into_iter().map(|l| IntervalsType::Window(l)));
    }
    if let Some(win_lens) = args.intervals_border {
        intervals_types.extend(win_lens.into_iter().map(|l| IntervalsType::Border(l)));
    }
    if let Some(seqs) = args.intervals_match {
        intervals_types.extend(seqs.into_iter().map(|mut s| {
            s.make_ascii_uppercase();
            IntervalsType::Match(s)
        }));
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    run(
        args.input,
        args.reference,
        args.stats_prefix,
        bin_types,
        intervals_types,
        args.name_column,
    );

    let duration = start_time.elapsed();
    println!("Run time (s): {}", duration.as_secs());
}

enum IntervalsType {
    Bed(Intervals),
    Homopolymer,
    Window(usize),
    Border(usize),
    Match(String),
}

#[derive(Parser)]
#[clap(author, version, about)]
struct Args {
    /// Input BAM file.
    input: String,

    /// Input reference FASTA file.
    reference: String,

    /// Prefix for output files that contain statistics.
    stats_prefix: String,

    /// Add column with a specific name in CSV outputs.
    #[clap(short, long)]
    name_column: Option<String>,

    /// Types of bins to use for per alignment stats.
    ///
    /// Each bin should be of the format <bin_type>:<step_size>.
    ///
    /// Supported bin types:
    /// q_len (read sequence length),
    /// subread_passes,
    /// mapq,
    /// mean_qual,
    /// gc_content,
    /// concordance_qv (phred scale Q-value)
    #[clap(short, long, min_values = 1)]
    bin_types: Option<Vec<String>>,

    /// Use intervals from a BED file.
    #[clap(long, min_values = 1)]
    intervals_bed: Option<Vec<String>>,

    /// Use homopolymer regions as intervals.
    #[clap(long)]
    intervals_hp: bool,

    /// Use fixed-width windows as intervals.
    #[clap(long, min_values = 1)]
    intervals_window: Option<Vec<usize>>,

    /// Use fixed-width window border regions as intervals.
    #[clap(long, min_values = 1)]
    intervals_border: Option<Vec<usize>>,

    /// Use regions that match any of the specified subsequences as intervals.
    #[clap(long, min_values = 1)]
    intervals_match: Option<Vec<String>>,

    /// Number of threads. Will be automatically determined if this is set to 0.
    #[clap(short, long, default_value_t = 0usize)]
    threads: usize,
}
