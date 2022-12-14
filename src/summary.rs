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

use std::fmt;

use fxhash::FxHashMap;

use ordered_float::OrderedFloat;

use crate::stats::*;

// important to ensure that summary stats are sorted so output order is deterministic

pub struct YieldSummary {
    name_column: Option<String>,
    /// (reads, bases)
    q_yield: [(usize, usize); 15],
}

impl YieldSummary {
    pub fn new(mut name_column: Option<String>) -> Self {
        if let Some(ref mut name) = name_column {
            name.push(',');
        }
        Self {
            name_column,
            q_yield: [(0usize, 0usize); 15],
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats) {
        if aln_stats.supplementary {
            return;
        }

        for i in 0..self.q_yield.len() {
            let min_q = i * 5;
            if aln_stats.concordance_qv >= (min_q as f64) {
                self.q_yield[i].0 += 1;
                self.q_yield[i].1 += aln_stats.q_len;
            }
        }
    }
}

impl fmt::Display for YieldSummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}min_empirical_q,yield_reads,yield_bases",
            if self.name_column.is_some() {
                "name,"
            } else {
                ""
            }
        )?;
        for i in 0..self.q_yield.len() {
            writeln!(
                f,
                "{}{},{},{}",
                self.name_column.as_ref().map(|n| n.as_str()).unwrap_or(""),
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
    name_column: Option<String>,
    pub total_alns: usize,
    matches: usize,
    mismatches: usize,
    non_hp_ins: usize,
    non_hp_del: usize,
    hp_ins: usize,
    hp_del: usize,
    gc_ins: usize,
    gc_del: usize,
    num_reads: usize,
}

impl IdentitySummary {
    pub fn new(mut name_column: Option<String>) -> Self {
        if let Some(ref mut name) = name_column {
            name.push(',');
        }
        Self {
            name_column,
            ..Default::default()
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats) {
        if aln_stats.supplementary {
            return;
        }

        self.matches += aln_stats.matches;
        self.mismatches += aln_stats.mismatches;
        self.non_hp_ins += aln_stats.non_hp_ins;
        self.non_hp_del += aln_stats.non_hp_del;
        self.hp_ins += aln_stats.hp_ins;
        self.hp_del += aln_stats.hp_del;
        self.gc_ins += aln_stats.gc_ins;
        self.gc_del += aln_stats.gc_del;
        self.num_reads += 1;
    }
}

impl fmt::Display for IdentitySummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}total_alns,primary_alns,identity,identity_qv,gap_compressed_identity,matches_per_kbp,mismatches_per_kbp,non_hp_ins_per_kbp,non_hp_del_per_kbp,hp_ins_per_kbp,hp_del_per_kbp", if self.name_column.is_some() { "name," } else { "" })?;
        let num_errors =
            self.mismatches + self.non_hp_ins + self.hp_ins + self.non_hp_del + self.hp_del;
        let num_bases = self.matches + self.mismatches + self.non_hp_del + self.hp_del;
        let id = (self.matches as f64) / ((self.matches + num_errors) as f64);
        let gc_id = (self.matches as f64)
            / ((self.matches + self.mismatches + self.gc_ins + self.gc_del) as f64);
        let per_kbp = |x| (x as f64) / (num_bases as f64) * 1000.0f64;
        writeln!(
            f,
            "{}{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
            self.name_column.as_ref().map(|n| n.as_str()).unwrap_or(""),
            self.total_alns,
            self.num_reads,
            id,
            concordance_qv(id, num_errors > 0),
            gc_id,
            per_kbp(self.matches),
            per_kbp(self.mismatches),
            per_kbp(self.non_hp_ins),
            per_kbp(self.non_hp_del),
            per_kbp(self.hp_ins),
            per_kbp(self.hp_del)
        )
    }
}

pub struct FeatureSummary {
    name_column: Option<String>,
    feature_stats: FxHashMap<String, FeatureStats>,
}

impl FeatureSummary {
    pub fn new(mut name_column: Option<String>) -> Self {
        if let Some(ref mut name) = name_column {
            name.push(',');
        }
        Self {
            name_column,
            feature_stats: FxHashMap::default(),
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats) {
        if aln_stats.supplementary {
            return;
        }

        for (&k, v) in &aln_stats.feature_stats {
            self.feature_stats
                .entry(k.to_owned())
                .or_insert_with(|| FeatureStats::default())
                .assign_add(v);
        }
    }
}

impl fmt::Display for FeatureSummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{}feature,intervals,identical_intervals,identity,identity_qv,mean_qual,bases_per_interval,matches_per_interval,mismatches_per_interval,non_hp_ins_per_interval,non_hp_del_per_interval,hp_ins_per_interval,hp_del_per_interval", if self.name_column.is_some() { "name," } else { "" })?;
        let mut v = self.feature_stats.iter().collect::<Vec<_>>();
        v.sort_by_key(|x| x.0);
        for (feature, stats) in v.into_iter() {
            let per_interval = |x| (x as f64) / (stats.overlaps as f64);
            let id = stats.identity();
            writeln!(
                f,
                "{}{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                self.name_column.as_ref().map(|n| n.as_str()).unwrap_or(""),
                feature.trim(),
                stats.overlaps,
                per_interval(stats.identical_overlaps),
                id,
                concordance_qv(id, id != 1.0),
                stats.mean_qual(),
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

pub struct CigarLenSummary {
    name_column: Option<String>,
    cigar_len_stats: FxHashMap<(usize, u8), usize>,
}

impl CigarLenSummary {
    pub fn new(mut name_column: Option<String>) -> Self {
        if let Some(ref mut name) = name_column {
            name.push(',');
        }
        Self {
            name_column,
            cigar_len_stats: FxHashMap::default(),
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats) {
        if aln_stats.supplementary {
            return;
        }

        for (&k, v) in &aln_stats.cigar_len_stats {
            *self.cigar_len_stats.entry(k).or_insert(0) += v;
        }
    }
}

impl fmt::Display for CigarLenSummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}cigar,length,count,length_count_per_cigar",
            if self.name_column.is_some() {
                "name,"
            } else {
                ""
            }
        )?;
        let mut v = self.cigar_len_stats.iter().collect::<Vec<_>>();
        v.sort_by_key(|(x, _)| (x.1, x.0));
        let mut total_cigars = [0usize; 128];
        for (cigar, &count) in &v {
            total_cigars[cigar.1 as usize] += count;
        }

        for (cigar, &count) in v.into_iter() {
            writeln!(
                f,
                "{}{},{},{},{:.6}",
                self.name_column.as_ref().map(|n| n.as_str()).unwrap_or(""),
                cigar.1 as char,
                cigar.0,
                count,
                (count as f64) / (total_cigars[cigar.1 as usize] as f64),
            )?;
        }
        Ok(())
    }
}

pub struct BinSummary {
    name_column: Option<String>,
    bin_maps: Vec<(BinType, FxHashMap<String, BinStats>)>,
}

impl BinSummary {
    pub fn new(mut name_column: Option<String>, bin_types: Vec<BinType>) -> Self {
        if let Some(ref mut name) = name_column {
            name.push(',');
        }
        let bin_maps = bin_types
            .into_iter()
            .map(|b| (b, FxHashMap::default()))
            .collect();
        Self {
            name_column,
            bin_maps,
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats) {
        if aln_stats.supplementary {
            return;
        }

        self.bin_maps.iter_mut().for_each(|(bin_type, bin_map)| {
            let bin = bin_type.get_bin(aln_stats);
            let bin_stats = BinStats::new(aln_stats);
            bin_map
                .entry(bin)
                .or_insert_with(|| BinStats::default())
                .assign_add(&bin_stats);
        });
    }
}

impl fmt::Display for BinSummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}bin_type,bin,num_reads,num_bases,identity,identity_qv,matches_per_kbp,mismatches_per_kbp,non_hp_ins_per_kbp,non_hp_del_per_kbp,hp_ins_per_kbp,hp_del_per_kbp",
            if self.name_column.is_some() {
                "name,"
            } else {
                ""
            }
        )?;
        for (bin_type, bin_map) in &self.bin_maps {
            let mut bins = bin_map.iter().collect::<Vec<_>>();
            bins.sort_by_key(|(b, _)| OrderedFloat(b.parse::<f64>().unwrap()));

            for (bin, stats) in bins {
                let per_kbp = |x| (x as f64) / (stats.num_bases() as f64) * 1000.0f64;
                let id = stats.identity();
                writeln!(
                    f,
                    "{}{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                    self.name_column.as_ref().map(|n| n.as_str()).unwrap_or(""),
                    bin_type,
                    bin,
                    stats.num_reads,
                    stats.num_bases(),
                    id,
                    concordance_qv(id, id != 1.0),
                    per_kbp(stats.matches),
                    per_kbp(stats.mismatches),
                    per_kbp(stats.non_hp_ins),
                    per_kbp(stats.non_hp_del),
                    per_kbp(stats.hp_ins),
                    per_kbp(stats.hp_del),
                )?;
            }
        }
        Ok(())
    }
}

pub struct QualScoreSummary {
    name_column: Option<String>,
    feature_qual: FxHashMap<String, QualScoreStats>,
}

impl QualScoreSummary {
    pub fn new(mut name_column: Option<String>) -> Self {
        if let Some(ref mut name) = name_column {
            name.push(',');
        }
        let mut feature_qual = FxHashMap::default();
        feature_qual.insert("all_alignments".to_owned(), QualScoreStats::default());
        Self {
            name_column,
            feature_qual,
        }
    }

    pub fn update(&mut self, aln_stats: &AlnStats) {
        if aln_stats.supplementary {
            return;
        }

        self.feature_qual
            .get_mut("all_alignments")
            .unwrap()
            .assign_add(&aln_stats.q_score_stats);

        for (&k, v) in &aln_stats.feature_stats {
            self.feature_qual
                .entry(k.to_owned())
                .or_insert_with(|| QualScoreStats::default())
                .assign_add(&v.q_score_stats);
        }
    }
}

impl fmt::Display for QualScoreSummary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "{}feature,qual_score,empirical_qv",
            if self.name_column.is_some() {
                "name,"
            } else {
                ""
            }
        )?;
        let mut v = self.feature_qual.iter().collect::<Vec<_>>();
        v.sort_by_key(|x| x.0);
        for (feature, stats) in v.into_iter() {
            for (i, qv) in stats.empirical_qv().into_iter() {
                writeln!(
                    f,
                    "{}{},{},{:.2}",
                    self.name_column.as_ref().map(|n| n.as_str()).unwrap_or(""),
                    feature.trim(),
                    i,
                    qv
                )?;
            }
        }
        Ok(())
    }
}
