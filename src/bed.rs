use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

use fxhash::{FxHashMap, FxHashSet};

use rust_lapper::{Interval, Lapper};

pub type FeatureInterval = Interval<usize, String>;

pub struct Intervals {
    intervals: FxHashMap<String, Lapper<usize, String>>,
    pub features: FxHashSet<String>,
}

impl Intervals {
    /// Create a new collection of intervals from a BED file.
    pub fn new(bed_path: &str) -> Self {
        let mut intervals = FxHashMap::default();
        let mut features = FxHashSet::default();
        let reader = BufReader::new(File::open(bed_path).expect("BED file not found."));

        for line in reader.lines() {
            let line = line.unwrap();
            let mut fields = line.split('\t');
            let chrom = fields.next().unwrap().to_owned();
            // convert to 1-indexed [start, stop)
            let start = fields.next().unwrap().parse::<usize>().unwrap() + 1;
            let stop = fields.next().unwrap().parse::<usize>().unwrap() + 1;
            let feature = fields.next().unwrap().to_owned();

            intervals
                .entry(chrom)
                .or_insert_with(|| Vec::<FeatureInterval>::new())
                .push(Interval {
                    start,
                    stop,
                    val: feature.clone(),
                });
            features.insert(feature);
        }

        Self {
            intervals: intervals
                .into_iter()
                .map(|(k, v)| (k, Lapper::new(v)))
                .collect(),
            features,
        }
    }

    /// Find intervals that intersect a given interval.
    pub fn find(&self, chrom: &str, start: usize, end: usize) -> Vec<&FeatureInterval> {
        self.intervals
            .get(chrom)
            .map(|x| x.find(start, end).collect())
            .unwrap_or_else(|| Vec::new())
    }
}
