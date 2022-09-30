use noodles::bam;
use noodles::core::Position;
use noodles::fasta;
use noodles::sam;

use sam::record::cigar::op::Kind;
use sam::record::data::field::Tag;

use fxhash::FxHashMap;

use std::fmt;
use std::str::FromStr;

use crate::bed::*;

/// Statistics for each alignment.
#[derive(Debug)]
pub struct AlnStats<'a> {
    pub read_name: String,
    pub q_len: usize,
    pub effective_cov: Option<f64>,
    pub subread_passes: Option<usize>,
    pub pred_concordance: Option<f64>,
    pub supplementary: bool,
    pub strand_rev: bool,
    pub mapq: u8,
    pub mean_qual: u8,
    pub read_len: usize,
    pub ref_cov: f64,
    pub gc_content: f64,
    pub concordance: f64,
    pub concordance_gc: f64,
    pub concordance_qv: f64,
    pub matches: usize,
    pub mismatches: usize,
    pub non_hp_ins: usize,
    pub non_hp_del: usize,
    pub hp_ins: usize,
    pub hp_del: usize,
    pub gc_ins: usize,
    pub gc_del: usize,
    pub feature_stats: FxHashMap<&'a str, FeatureStats>,
    pub cigar_len_stats: FxHashMap<(usize, u8), usize>,
}

/// Per-read attributes that can be binned.
#[derive(Copy, Clone)]
pub enum BinType {
    QLen(usize),
    SubreadPasses(usize),
    MapQ(u8),
    MeanQual(u8),
    GcContent(f64),
    ConcordanceQv(f64),
}

impl BinType {
    pub fn get_bin(&self, a: &AlnStats) -> String {
        match self {
            Self::QLen(step) => format!("{}", a.q_len / step * step),
            Self::SubreadPasses(step) => format!(
                "{}",
                a.subread_passes.expect("Subread passes not found!") / step * step
            ),
            Self::MapQ(step) => format!("{}", a.mapq / step * step),
            Self::MeanQual(step) => format!("{}", a.mean_qual / step * step),
            Self::GcContent(step) => format!("{:.6}", (a.gc_content / step).floor() * step),
            Self::ConcordanceQv(step) => format!("{:.2}", (a.concordance_qv / step).floor() * step),
        }
    }
}

impl fmt::Display for BinType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::QLen(step) => write!(f, "q_len:{}", step),
            Self::SubreadPasses(step) => write!(f, "subread_passes:{}", step),
            Self::MapQ(step) => write!(f, "mapq:{}", step),
            Self::MeanQual(step) => write!(f, "mean_qual:{}", step),
            Self::GcContent(step) => write!(f, "gc_content:{}", step),
            Self::ConcordanceQv(step) => write!(f, "concordance_qv:{}", step),
        }
    }
}

impl FromStr for BinType {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut split = s.split(':');
        let a = split
            .next()
            .expect("Bin type not found! Expected <bin_type>:<step_size>");
        let b = split
            .next()
            .expect("Step size not found! Expected <bin_type>:<step_size>");

        use BinType::*;
        match a {
            "q_len" => Ok(QLen(b.parse::<_>().unwrap())),
            "subread_passes" => Ok(SubreadPasses(b.parse::<_>().unwrap())),
            "mapq" => Ok(MapQ(b.parse::<_>().unwrap())),
            "mean_qual" => Ok(MeanQual(b.parse::<_>().unwrap())),
            "gc_content" => Ok(GcContent(b.parse::<_>().unwrap())),
            "concordance_qv" => Ok(ConcordanceQv(b.parse::<_>().unwrap())),
            _ => Err("Invalid stat to bin across!".into()),
        }
    }
}

/// Statistics for each bin.
#[derive(Debug, Default)]
pub struct BinStats {
    pub num_reads: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub non_hp_ins: usize,
    pub non_hp_del: usize,
    pub hp_ins: usize,
    pub hp_del: usize,
}

impl BinStats {
    pub fn new(stats: &AlnStats) -> Self {
        Self {
            num_reads: 1,
            matches: stats.matches,
            mismatches: stats.mismatches,
            non_hp_ins: stats.non_hp_ins,
            non_hp_del: stats.non_hp_del,
            hp_ins: stats.hp_ins,
            hp_del: stats.hp_del,
        }
    }

    pub fn assign_add(&mut self, o: &Self) {
        self.num_reads += o.num_reads;
        self.matches += o.matches;
        self.mismatches += o.mismatches;
        self.non_hp_ins += o.non_hp_ins;
        self.non_hp_del += o.non_hp_del;
        self.hp_ins += o.hp_ins;
        self.hp_del += o.hp_del;
    }

    pub fn num_bases(&self) -> usize {
        self.matches + self.mismatches + self.non_hp_del + self.hp_del
    }

    pub fn num_errors(&self) -> usize {
        self.mismatches + self.non_hp_ins + self.hp_ins + self.non_hp_del + self.hp_del
    }

    pub fn identity(&self) -> f64 {
        (self.matches as f64) / ((self.matches + self.num_errors()) as f64)
    }
}

/// Statistics for each interval feature.
#[derive(Debug, Default)]
pub struct FeatureStats {
    pub overlaps: usize,
    pub identical_overlaps: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub non_hp_ins: usize,
    pub non_hp_del: usize,
    pub hp_ins: usize,
    pub hp_del: usize,
}

impl FeatureStats {
    pub fn assign_add(&mut self, o: &Self) {
        self.overlaps += o.overlaps;
        self.identical_overlaps += o.identical_overlaps;
        self.matches += o.matches;
        self.mismatches += o.mismatches;
        self.non_hp_ins += o.non_hp_ins;
        self.non_hp_del += o.non_hp_del;
        self.hp_ins += o.hp_ins;
        self.hp_del += o.hp_del;
    }

    pub fn num_bases(&self) -> usize {
        self.matches + self.mismatches + self.non_hp_del + self.hp_del
    }

    pub fn num_errors(&self) -> usize {
        self.mismatches + self.non_hp_ins + self.hp_ins + self.non_hp_del + self.hp_del
    }

    pub fn identity(&self) -> f64 {
        (self.matches as f64) / ((self.matches + self.num_errors()) as f64)
    }
}

impl<'a> AlnStats<'a> {
    pub fn from_record(
        references: &sam::header::ReferenceSequences,
        reference_seqs: &FxHashMap<String, fasta::Record>,
        r: &bam::lazy::Record,
        intervals: &[&'a FeatureInterval],
    ) -> Self {
        // note: avoid copying data (especially sequence/quality scores) since they are large
        let sequence = sam::record::Sequence::try_from(r.sequence()).unwrap();
        let q_scores = sam::record::QualityScores::try_from(r.quality_scores()).unwrap();
        let flags = r.flags().unwrap();
        let data = sam::record::Data::try_from(r.data()).unwrap();
        let ec_tag = Tag::try_from(*b"ec").unwrap();
        let ec = data
            .get(ec_tag)
            .map(|f| f.value().as_float().unwrap() as f64);
        let np_tag: Tag = Tag::try_from(*b"np").unwrap();
        let np = data
            .get(np_tag)
            .map(|f| f.value().as_int().unwrap() as usize);
        let rq_tag: Tag = Tag::try_from(*b"rq").unwrap();
        let rq = data
            .get(rq_tag)
            .map(|f| f.value().as_float().unwrap() as f64);

        let mut res = AlnStats {
            read_name: r.read_name().unwrap().unwrap().to_string(),
            q_len: sequence.len(),
            effective_cov: ec,
            subread_passes: np,
            pred_concordance: rq,
            supplementary: flags.is_supplementary(),
            strand_rev: flags.is_reverse_complemented(),
            mapq: u8::from(r.mapping_quality().unwrap().unwrap()),
            mean_qual: mean_qual(q_scores.as_ref()),
            // fill in the rest afterwards
            read_len: 0,
            ref_cov: 0.0,
            gc_content: 0.0,
            concordance: 0.0,
            concordance_gc: 0.0,
            concordance_qv: 0.0,
            matches: 0,
            mismatches: 0,
            non_hp_ins: 0,
            non_hp_del: 0,
            hp_ins: 0,
            hp_del: 0,
            gc_ins: 0,
            gc_del: 0,
            feature_stats: FxHashMap::default(),
            cigar_len_stats: FxHashMap::default(),
        };

        let mut interval_has_error = vec![false; intervals.len()];
        for i in intervals {
            res.feature_stats
                .entry(&i.val)
                .or_insert_with(|| FeatureStats::default())
                .overlaps += 1;
        }

        let mut ref_pos = usize::from(r.alignment_start().unwrap().unwrap());
        let mut query_pos = 1;
        let mut interval_start_idx = 0;
        let curr_ref_name = references[r.reference_sequence_id().unwrap().unwrap()]
            .name()
            .to_string();
        let curr_ref_seq = reference_seqs[&curr_ref_name].sequence();
        let mut curr_features = Vec::new();
        let mut curr_interval_idxs = Vec::new();

        let mut intervals_have_error = |v: &[usize]| {
            v.iter()
                .for_each(|&interval_idx| interval_has_error[interval_idx] = true);
        };

        // count mismatches, indels, and homopolymers
        let cigar = sam::record::Cigar::try_from(r.cigar()).unwrap();
        for op in cigar.iter() {
            for _i in 0..op.len() {
                // skip intervals that cannot overlap the current reference position
                while interval_start_idx < intervals.len()
                    && ref_pos >= intervals[interval_start_idx].stop
                {
                    interval_start_idx += 1;
                }
                // find the intervals that overlap the current reference position
                let mut interval_idx = interval_start_idx;
                curr_features.clear();
                curr_interval_idxs.clear();
                while interval_idx < intervals.len() && ref_pos >= intervals[interval_idx].start {
                    if ref_pos < intervals[interval_idx].stop {
                        // get feature names of the overlapping intervals
                        curr_features.push(intervals[interval_idx].val.as_str());
                        curr_interval_idxs.push(interval_idx);
                    }
                    interval_idx += 1;
                }

                match op.kind() {
                    Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match => {
                        let c = curr_ref_seq[Position::new(ref_pos).unwrap()].to_ascii_uppercase();
                        let is_match = op.kind() == Kind::SequenceMatch
                            || (op.kind() == Kind::Match
                                && c == u8::from(sequence[Position::new(query_pos).unwrap()])
                                    .to_ascii_uppercase());
                        if is_match {
                            res.matches += 1;
                            curr_features
                                .iter()
                                .for_each(|f| res.feature_stats.get_mut(f).unwrap().matches += 1);
                        } else {
                            res.mismatches += 1;
                            curr_features.iter().for_each(|f| {
                                res.feature_stats.get_mut(f).unwrap().mismatches += 1
                            });
                            intervals_have_error(&curr_interval_idxs);
                        }
                        if c == b'C' || c == b'G' {
                            res.gc_content += 1.0;
                        }
                        query_pos += 1;
                        ref_pos += 1;
                    }
                    Kind::Insertion => {
                        // can be computed without looping through the number of insertions
                        // this does not modify ref_pos
                        let before_ins =
                            curr_ref_seq[Position::new(ref_pos).unwrap()].to_ascii_uppercase();
                        let after_ins = curr_ref_seq
                            .get(Position::new(ref_pos + 1).unwrap())
                            .unwrap_or(&b'?')
                            .to_ascii_uppercase();
                        let query_ins = &sequence[Position::new(query_pos).unwrap()
                            ..Position::new(query_pos + op.len()).unwrap()];
                        let hp_before = query_ins
                            .iter()
                            .map(|&c| u8::from(c).to_ascii_uppercase())
                            .all(|c| c == before_ins);
                        let hp_after = query_ins
                            .iter()
                            .map(|&c| u8::from(c).to_ascii_uppercase())
                            .all(|c| c == after_ins);
                        if hp_before || hp_after {
                            res.hp_ins += op.len();
                            curr_features
                                .iter()
                                .for_each(|f| res.feature_stats.get_mut(f).unwrap().hp_ins += 1);
                        } else {
                            res.non_hp_ins += op.len();
                            curr_features.iter().for_each(|f| {
                                res.feature_stats.get_mut(f).unwrap().non_hp_ins += 1
                            });
                        }
                        intervals_have_error(&curr_interval_idxs);
                        query_pos += op.len();
                        break;
                    }
                    Kind::Deletion => {
                        let before_curr = curr_ref_seq
                            .get(Position::new(ref_pos - 1).unwrap())
                            .unwrap_or(&b'?')
                            .to_ascii_uppercase();
                        let after_curr = curr_ref_seq
                            .get(Position::new(ref_pos + 1).unwrap())
                            .unwrap_or(&b'?')
                            .to_ascii_uppercase();
                        let curr =
                            curr_ref_seq[Position::new(ref_pos).unwrap()].to_ascii_uppercase();
                        if curr == b'C' || curr == b'G' {
                            res.gc_content += 1.0;
                        }
                        let hp = curr == before_curr || curr == after_curr;
                        if hp {
                            res.hp_del += 1;
                            curr_features
                                .iter()
                                .for_each(|f| res.feature_stats.get_mut(f).unwrap().hp_del += 1);
                        } else {
                            res.non_hp_del += 1;
                            curr_features.iter().for_each(|f| {
                                res.feature_stats.get_mut(f).unwrap().non_hp_del += 1
                            });
                        }
                        intervals_have_error(&curr_interval_idxs);
                        ref_pos += 1;
                    }
                    Kind::SoftClip => {
                        // does not require looping through the number of soft clips
                        query_pos += op.len();
                        break;
                    }
                    Kind::HardClip => {
                        // does not require looping through the number of hard clips
                        break;
                    }
                    Kind::Skip => {
                        // does not require looping through the number of skip operations
                        ref_pos += op.len();
                        break;
                    }
                    _ => panic!("Unexpected CIGAR operation: {}", op),
                }
            }

            // gap compressed
            match op.kind() {
                Kind::SequenceMatch => {
                    *res.cigar_len_stats.entry((op.len(), b'=')).or_insert(0) += 1;
                }
                Kind::SequenceMismatch => {
                    *res.cigar_len_stats.entry((op.len(), b'X')).or_insert(0) += 1;
                }
                Kind::Match => {
                    *res.cigar_len_stats.entry((op.len(), b'M')).or_insert(0) += 1;
                }
                Kind::Insertion => {
                    *res.cigar_len_stats.entry((op.len(), b'I')).or_insert(0) += 1;
                    res.gc_ins += 1;
                }
                Kind::Deletion => {
                    *res.cigar_len_stats.entry((op.len(), b'D')).or_insert(0) += 1;
                    res.gc_del += 1;
                }
                _ => (),
            }
        }

        let errors = res.mismatches + res.non_hp_ins + res.non_hp_del + res.hp_ins + res.hp_del;
        res.read_len = res.matches + res.mismatches + res.non_hp_del + res.hp_del;
        res.ref_cov = (res.read_len as f64) / (curr_ref_seq.len() as f64);
        res.gc_content /= res.read_len as f64;
        res.concordance = (res.matches as f64) / ((res.matches + errors) as f64);
        res.concordance_gc = (res.matches as f64)
            / ((res.matches + res.mismatches + res.gc_ins + res.gc_del) as f64);
        res.concordance_qv = concordance_qv(res.concordance, errors > 0);

        for (i, &has_error) in intervals.iter().zip(&interval_has_error) {
            if !has_error {
                res.feature_stats
                    .get_mut(i.val.as_str())
                    .unwrap()
                    .identical_overlaps += 1;
            }
        }

        res
    }

    pub fn header() -> &'static str {
        "read,read_length,effective_coverage,subread_passes,predicted_concordance,alignment_type,strand,alignment_mapq,mean_quality,aligned_read_length,reference_coverage,gc_content,concordance,gap_compressed_concordance,concordance_qv,mismatches,non_hp_ins,non_hp_del,hp_ins,hp_del"
    }

    pub fn to_csv(&self) -> String {
        let supp_str = if self.supplementary {
            "supplementary"
        } else {
            "primary"
        };
        let strand_str = if self.strand_rev { "-" } else { "+" };
        let ec = self
            .effective_cov
            .map(|x| format!("{:.2}", x))
            .unwrap_or_else(|| String::new());
        let np = self
            .subread_passes
            .map(|x| format!("{}", x))
            .unwrap_or_else(|| String::new());
        let rq = self
            .pred_concordance
            .map(|x| format!("{:.6}", x))
            .unwrap_or_else(|| String::new());
        format!(
            "{},{},{:.2},{},{:.6},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.2},{},{},{},{},{}",
            self.read_name,
            self.q_len,
            ec,
            np,
            rq,
            supp_str,
            strand_str,
            self.mapq,
            self.mean_qual,
            self.read_len,
            self.ref_cov,
            self.gc_content,
            self.concordance,
            self.concordance_gc,
            self.concordance_qv,
            self.mismatches,
            self.non_hp_ins,
            self.non_hp_del,
            self.hp_ins,
            self.hp_del
        )
    }
}

/// Compute the Phred scale Q-value for a certain concordance/identity.
///
/// Perfect match will have a Q-value of 75.
pub fn concordance_qv(concordance: f64, has_errors: bool) -> f64 {
    if has_errors {
        -10.0f64 * (1.0f64 - concordance).log10()
    } else {
        75.0f64
    }
}

fn mean_qual(q_scores: &[sam::record::quality_scores::Score]) -> u8 {
    let sum_q = q_scores
        .iter()
        .map(|&q| 10.0f64.powf(-(u8::from(q) as f64) / 10.0f64))
        .sum::<f64>();
    (-10.0f64 * (sum_q / (q_scores.len() as f64)).log10()).round() as u8
}
