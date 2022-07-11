use std::fmt;

use crate::stats::*;

pub struct YieldSummary {
    prefix: String,
    /// (reads, bases)
    q_yield: [(usize, usize); 10],
}

impl YieldSummary {
    pub fn new(prefix: String) -> Self {
        Self {
            prefix,
            q_yield: [(0usize, 0usize); 10],
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats) {
        if aln_stats.supplementary {
            return;
        }

        for i in 0..self.q_yield.len() {
            let min_q = i * 5;
            if aln_stats.concordance_qv >= (min_q as f32) {
                self.q_yield[i].0 += 1;
                self.q_yield[i].1 += aln_stats.q_len;
            }
        }
    }
}

impl fmt::Display for YieldSummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "prefix,min_empirical_Q,yield_reads,yield_bases")?;
        for i in 0..self.q_yield.len() {
            writeln!(
                f,
                "{},{},{},{}",
                self.prefix,
                i * 5,
                self.q_yield[i].0,
                self.q_yield[i].1
            )?;
        }
        Ok(())
    }
}

#[derive(Default)]
pub struct IdentitySummary {
    prefix: String,
    matches: usize,
    mismatches: usize,
    ins: usize,
    del: usize,
    gc_ins: usize,
    gc_del: usize,
}

impl IdentitySummary {
    pub fn new(prefix: String) -> Self {
        Self {
            prefix,
            ..Default::default()
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats) {
        if aln_stats.supplementary {
            return;
        }

        self.matches += aln_stats.matches;
        self.mismatches += aln_stats.mismatches;
        self.ins += aln_stats.non_hp_ins + aln_stats.hp_ins;
        self.del += aln_stats.non_hp_del + aln_stats.hp_del;
        self.gc_ins += aln_stats.gc_ins;
        self.gc_del += aln_stats.gc_del;
    }
}

impl fmt::Display for IdentitySummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "prefix,identity,gap_compressed_identity")?;
        let id =
            (self.matches as f32) / ((self.matches + self.mismatches + self.ins + self.del) as f32);
        let gc_id = (self.matches as f32)
            / ((self.matches + self.mismatches + self.gc_ins + self.gc_del) as f32);
        writeln!(f, "{},{},{}", self.prefix, id, gc_id)
    }
}
