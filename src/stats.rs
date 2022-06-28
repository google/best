use noodles_sam as sam;
use noodles_fasta as fasta;

#[derive(Debug)]
struct AlnStats {
    read_name: String,
    q_len: usize,
    effective_cov: f64,
    subread_passes: usize,
    pred_concordance: f64,
    supplementary: bool,
    mapq: usize,
    interval_q_len: usize,
    concordance: f64,
    concordance_qv: f64,
    mismatches: usize,
    non_hp_ins: usize,
    non_hp_del: usize,
    hp_ins: usize,
    hp_del: usize,
}

impl AlnStats {
    fn from_record(references: &fasta::Repository, r: &sam::Record) -> Self {
        let mut res = AlnStats {
            read_name: r.read_name().unwrap().to_string(),
            q_len: r.sequence().len(),
            effective_cov: r.data().get("ec".into()).unwrap_or(0),
            subread_passes: r.data().get("np".into()).unwrap_or(0),
            pred_concordance: r.data().get("rq".into()).unwrap_or(0),
            supplementary: r.flags().is_supplementary(),
            mapq: r.mapping_quality().unwrap(),
            interval_q_len: r.sequence().len(),
            concordance: 0,
            concordance_qv: 0,
            mismatches: 0,
            non_hp_ins: 0,
            non_hp_del: 0,
            hp_ins: 0,
            hp_del: 0,
        };

        let mut matches = 0;

        for op in cigar {
            match op.kind {
                Kind::SequenceMatch => matches += 1,
                Kind::SequenceMismatch => args.mismatches += 1,
                Kind::Insertion => {

                },
                Kind::Deletion => {

                },
                _ => panic!("Unexpected CIGAR operation!"),
            }
        }

        let errors = args.mismatches + args.non_hp_ins + args.non_hp_del + args.hp_ins + args.hp_del;
        res.concordance = (matches as f64) / ((matches + errors) as f64);
        res.concordance_qv = concordance_qv(res.concordance, res.q_len, errors > 0);

        res
    }

    fn header() -> &str {
        "#read,readLengthBp,effectiveCoverage,subreadPasses,predictedConcordance,alignmentType,alignmentMapq,hcReadLengthBp,concordance,concordanceQv,mismatchBp,nonHpInsertionBp,nonHpDeletionBp,hpInsertionBp,hpDeletionBp"
    }

    fn to_csv(&self) -> String {
        let supp_str = if self.supplementary { "Supplementary" } else { "Primary" };
        format!("{},{},{:.2f},{},{:.6f},{},{},{},{:.6f},{:.2f},{},{},{},{},{}", self.read_name, self.q_len, self.effective_cov, self.subread_passes, self.pred_concordance, supp_str, self.mapq, self.interval_q_len, self.concordance, self.concordance_qv, self.mismatches, self.non_hp_ins, self.non_hp_del, self.hp_ins, self.hp_del)
    }
}

fn concordance_qv(concordance: f64, read_len: usize, has_errors: bool) -> f64 {
    let qv_cap = 10.0f64 * ((read_len as f64) + 1.0f64).log10();
    if has_errors {
        qv_cap.min(-10.0f64 * (1.0f64 - concordance).log10())
    } else {
        qv_cap
    }
}
