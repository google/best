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
            let feature = fields.next().unwrap_or("none").to_owned();

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
