use std::fmt;

use fxhash::FxHashMap;

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
    hp_ins: usize,
    hp_del: usize,
    gc_ins: usize,
    gc_del: usize,
    num_reads: usize,
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
        self.hp_ins += aln_stats.hp_ins;
        self.hp_del += aln_stats.hp_del;
        self.gc_ins += aln_stats.gc_ins;
        self.gc_del += aln_stats.gc_del;
        self.num_reads += 1;
    }
}

impl fmt::Display for IdentitySummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "prefix,identity,gap_compressed_identity,matches_per_read,mismatches_per_read,ins_per_read,del_per_read,hp_ins_per_read,hp_del_per_read")?;
        let id =
            (self.matches as f32) / ((self.matches + self.mismatches + self.ins + self.del) as f32);
        let gc_id = (self.matches as f32)
            / ((self.matches + self.mismatches + self.gc_ins + self.gc_del) as f32);
        let per_read = |x| (x as f64) / (self.num_reads as f64);
        writeln!(
            f,
            "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
            self.prefix,
            id,
            gc_id,
            per_read(self.matches),
            per_read(self.mismatches),
            per_read(self.ins),
            per_read(self.del),
            per_read(self.hp_ins),
            per_read(self.hp_del)
        )
    }
}

pub struct FeatureSummary<'a> {
    prefix: String,
    feature_stats: FxHashMap<&'a str, FeatureStats>,
}

impl<'a> FeatureSummary<'a> {
    pub fn new(prefix: String) -> Self {
        Self {
            prefix,
            feature_stats: FxHashMap::default(),
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats<'a>) {
        if aln_stats.supplementary {
            return;
        }

        for (k, v) in &aln_stats.feature_stats {
            self.feature_stats
                .entry(k)
                .or_insert_with(|| FeatureStats::default())
                .assign_add(v);
        }
    }
}

impl<'a> fmt::Display for FeatureSummary<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "prefix,feature,bases_per_interval,matches_per_interval,mismatches_per_interval,non_hp_ins_per_interval,non_hp_del_per_interval,hp_ins_per_interval,hp_del_per_interval")?;
        let mut v = self.feature_stats.iter().collect::<Vec<_>>();
        v.sort_by_key(|x| x.0);
        for (feature, stats) in v.into_iter() {
            let per_interval = |x| (x as f64) / (stats.overlaps as f64);
            writeln!(
                f,
                "{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                self.prefix,
                feature,
                per_interval(stats.num_bases()),
                per_interval(stats.matches),
                per_interval(stats.mismatches),
                per_interval(stats.non_hp_ins),
                per_interval(stats.non_hp_del),
                per_interval(stats.hp_ins),
                per_interval(stats.hp_del)
            )?;
        }
        Ok(())
    }
}
