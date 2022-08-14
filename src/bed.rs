use fxhash::{FxHashMap, FxHashSet};

use rust_lapper::{Lapper, Interval};

pub struct Intervals {
    intervals: FxHashMap<String, Lapper<usize, String>>,
    pub features: FxHashSet<String>,
}

impl Intervals {
    pub fn new(bed_path: &str) -> Self {
        let mut intervals = FxHashMap::default();
        let mut features = FxHashSet::default();
        let mut reader = BufReader::new(File::open(bed_path).expect("BED file not found."));

        for line in reader.lines() {
            let fields = line.split('\t');
            let chrom = line.next().unwrap().to_owned();
            let start = line.next().unwrap().parse::<usize>();
            let stop = line.next().unwrap().parse::<usize>();
            let feature = line.next().unwrap().to_owned();

            intervals.entry(&chrom).or_insert_with(|| Vec<Interval>::new()).push(Interval {
                start,
                stop,
                val: feature,
            });
            features.insert(feature);
        }

        Self {
            intervals: intervals.into_iter().map(|k, v| (k, Lapper::new(v))).collect(),
            features,
        }
    }

    pub fn find(&self, chrom: &str, start: usize, end: usize) -> Vec<&Interval> {
        self.intervals.get(chrom).map(|x| x.find(start, end).collect()).unwrap_or_else(|| Vec::new())
    }
}
