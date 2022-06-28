use clap::Parser;

use rayon;

use noodles_bam as bam;
use noodles_fasta as fasta;

mod stats;
use stats::*;

fn run(input_path: String, reference_path: String, aln_stats_path: String) {
    let ref_index = fasta::index(reference_path);
    let ref_reader = fasta::Reader::new(BufReader::new(File::open(reference_path).unwrap()));
    let adapter = fasta::repository::adapters::IndexedReader::new(ref_reader, ref_index);
    let references = fasta::Repository::new(adapter);

    let mut reader = bam::Reader::new(File::open(input_path).unwrap());
    reader.read_header().unwrap();
    reader.read_reference_sequences().unwrap();

    let mut aln_stats_writer = File::create(aln_stats_path).unwrap();
    write!(aln_stats_writer, AlnStats::header());
    let aln_stats_writer = Mutex::new(aln_stats_writer);

    reader
        .records()
        .into_par_iter()
        .map(|r| r.unwrap())
        .filter(|r| !r.flags().is_unmapped() && !r.flags().is_secondary())
        .for_each(|record| {
            let stats = AlnStats::from_record(&references, record);

            {
                let writer = aln_stats_writer.lock();
                write!(writer, "{}\n", stats.to_csv()).unwrap();
            }
        });
}

fn main() {
    let args = Args::parse();

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    run(args.input, args.reference, args.aln_stats);
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
