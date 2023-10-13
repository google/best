// Copyright (c) 2022 Google LLC
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
// the Software, and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
// IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

use clap::Parser;

use rayon::prelude::*;

use noodles::{bam, fasta, sam};

use fxhash::FxHashMap;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
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

const PER_ALN_STATS_NAME: &str = "per_aln_stats.csv.gz";
const YIELD_STATS_NAME: &str = "summary_yield_stats.csv";
const IDENTITY_STATS_NAME: &str = "summary_identity_stats.csv";
const FEATURE_STATS_NAME: &str = "summary_feature_stats.csv";
const CIGAR_STATS_NAME: &str = "summary_cigar_stats.csv";
const BIN_STATS_NAME: &str = "summary_bin_stats.csv";
const QUAL_SCORE_STATS_NAME: &str = "summary_qual_score_stats.csv";

fn run(
    input_path: String,
    reference_path: String,
    stats_prefix: String,
    bin_types: Option<Vec<BinType>>,
    intervals_types: Vec<IntervalsType>,
    name_column: Option<String>,
    output_per_aln_stats: bool,
) {
    // read reference sequences from fasta file
    let mut ref_reader = {
        let f = File::open(&reference_path).unwrap();
        let r: Box<dyn Read> = if reference_path.ends_with(".gz") {
            Box::new(GzDecoder::new(f))
        } else {
            Box::new(f)
        };
        fasta::Reader::new(BufReader::new(r))
    };
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
    let aln_stats_writer = if output_per_aln_stats {
        let mut w = GzEncoder::new(
            BufWriter::new(File::create(&aln_stats_path).unwrap()),
            flate2::Compression::default(),
        );
        write!(
            w,
            "{}{}\n",
            if name_column.is_some() { "name," } else { "" },
            AlnStats::header()
        )
        .unwrap();
        Some(Mutex::new(w))
    } else {
        None
    };

    let summary_yield = Mutex::new(YieldSummary::new(name_column.clone()));
    let summary_identity = Mutex::new(IdentitySummary::new(name_column.clone()));
    let summary_features = if intervals_types.is_empty() {
        None
    } else {
        Some(Mutex::new(FeatureSummary::new(name_column.clone())))
    };
    let summary_cigars = Mutex::new(CigarLenSummary::new(name_column.clone()));
    let summary_bins = bin_types.map(|b| Mutex::new(BinSummary::new(name_column.clone(), b)));
    let summary_qual_score = Mutex::new(QualScoreSummary::new(name_column.clone()));
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

            let strand_rev = flags.is_reverse_complemented();
            let aln_ref = references[record.reference_sequence_id().unwrap().unwrap()]
                .name()
                .as_str();
            if !reference_seqs.contains_key(aln_ref) {
                panic!(
                    "{} is not found in the input reference sequence names!",
                    aln_ref
                );
            }
            // convert to one-indexed [aln_start, aln_end)
            let aln_start = usize::from(record.alignment_start().unwrap().unwrap());
            let aln_end = aln_start
                + sam::record::Cigar::try_from(record.cigar())
                    .unwrap()
                    .alignment_span();
            // get all the intervals relevant for the current alignment record
            let mut intervals_vec = Vec::new();
            let mut overlap_intervals = Vec::new();
            intervals_types
                .iter()
                .for_each(|intervals_type| match intervals_type {
                    IntervalsType::Homopolymer => intervals_vec.extend(find_homopolymers(
                        reference_seqs[aln_ref].sequence(),
                        aln_start,
                        aln_end,
                        strand_rev,
                    )),
                    IntervalsType::Window(win_len) => intervals_vec
                        .extend(get_windows(aln_start, aln_end, *win_len, false, strand_rev)),
                    IntervalsType::WindowPos(win_len) => intervals_vec
                        .extend(get_windows(aln_start, aln_end, *win_len, true, strand_rev)),
                    IntervalsType::Border(win_len) => {
                        intervals_vec.extend(get_borders(aln_start, aln_end, *win_len, strand_rev))
                    }
                    IntervalsType::Match(seq) => intervals_vec.extend(get_matches(
                        reference_seqs[aln_ref].sequence(),
                        aln_start,
                        aln_end,
                        seq,
                        strand_rev,
                    )),
                    IntervalsType::Bed(intervals) => {
                        overlap_intervals.extend(intervals.find(aln_ref, aln_start, aln_end))
                    }
                });
            overlap_intervals.extend(&intervals_vec);
            overlap_intervals.sort();

            let stats =
                AlnStats::from_record(&references, &reference_seqs, &record, &overlap_intervals);

            summary_yield.lock().unwrap().update(&stats);
            summary_identity.lock().unwrap().update(&stats);
            summary_features
                .as_ref()
                .map(|f| f.lock().unwrap().update(&stats));
            summary_cigars.lock().unwrap().update(&stats);
            summary_bins
                .as_ref()
                .map(|b| b.lock().unwrap().update(&stats));
            summary_qual_score.lock().unwrap().update(&stats);

            if let Some(ref w) = aln_stats_writer {
                let mut w = w.lock().unwrap();
                if let Some(ref name) = name_column {
                    write!(w, "{},{}\n", name, stats.to_csv()).unwrap();
                } else {
                    write!(w, "{}\n", stats.to_csv()).unwrap();
                }
            }
        });

    write_summary(
        summary_yield.into_inner().unwrap(),
        &stats_prefix,
        YIELD_STATS_NAME,
    );

    summary_identity.lock().unwrap().total_alns = total_alns.into_inner();
    write_summary(
        summary_identity.into_inner().unwrap(),
        &stats_prefix,
        IDENTITY_STATS_NAME,
    );

    if let Some(f) = summary_features {
        write_summary(f.into_inner().unwrap(), &stats_prefix, FEATURE_STATS_NAME);
    }

    write_summary(
        summary_cigars.into_inner().unwrap(),
        &stats_prefix,
        CIGAR_STATS_NAME,
    );

    if let Some(b) = summary_bins {
        write_summary(b.into_inner().unwrap(), &stats_prefix, BIN_STATS_NAME);
    }

    write_summary(
        summary_qual_score.into_inner().unwrap(),
        &stats_prefix,
        QUAL_SCORE_STATS_NAME,
    );
}

fn write_summary<D: std::fmt::Display>(s: D, prefix: &str, name: &str) {
    let summary_path = format!("{}.{}", prefix, name);
    let mut summary_writer = File::create(&summary_path).unwrap();
    write!(summary_writer, "{}", s).unwrap();
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
    if let Some(win_lens) = args.intervals_window_pos {
        intervals_types.extend(win_lens.into_iter().map(|l| IntervalsType::WindowPos(l)));
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
        !args.no_per_aln_stats,
    );

    let duration = start_time.elapsed();
    println!("Run time (s): {}", duration.as_secs());
}

enum IntervalsType {
    Bed(Intervals),
    Homopolymer,
    Window(usize),
    WindowPos(usize),
    Border(usize),
    Match(String),
}

#[derive(Parser)]
#[clap(author, version, about)]
struct Args {
    /// Input BAM file.
    input: String,

    /// Input reference FASTA file. Can be gzipped.
    reference: String,

    /// Prefix for output files that contain statistics.
    stats_prefix: String,

    /// Add column with a specific name in CSV outputs.
    #[clap(short, long)]
    name_column: Option<String>,

    /// Turn off outputting per alignment stats.
    #[clap(long)]
    no_per_aln_stats: bool,

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
    ///
    /// The BED file should have the columns chrom, start, stop, and feature.
    /// The feature column is optional.
    ///
    /// This allows stats to be gathered separately for different types of intervals.
    /// Note that all intervals are on the reference, not the reads.
    #[clap(long, min_values = 1)]
    intervals_bed: Option<Vec<String>>,

    /// Use homopolymer regions in the reference as intervals.
    #[clap(long)]
    intervals_hp: bool,

    /// Use fixed-width nonoverlapping windows as intervals.
    ///
    /// This is used to specify the window widths.
    #[clap(long, min_values = 1)]
    intervals_window: Option<Vec<usize>>,

    /// Use fixed-width nonoverlapping windows with positions as intervals.
    ///
    /// This is used to specify the window widths.
    #[clap(long, min_values = 1)]
    intervals_window_pos: Option<Vec<usize>>,

    /// Use fixed-width nonoverlapping window border regions as intervals.
    ///
    /// This is used to specify the window widths.
    #[clap(long, min_values = 1)]
    intervals_border: Option<Vec<usize>>,

    /// Use regions in the reference that match any of the specified subsequences as intervals.
    #[clap(long, min_values = 1)]
    intervals_match: Option<Vec<String>>,

    /// Number of threads. Will be automatically determined if this is set to 0.
    #[clap(short, long, default_value_t = 0usize)]
    threads: usize,
}
