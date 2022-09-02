use noodles::core::Position;
use noodles::fasta;

use crate::bed::FeatureInterval;

pub fn find_homopolymers(
    seq: &fasta::record::Sequence,
    start: usize,
    end: usize,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();
    let mut hp_len = 0;
    let mut prev = b'?';

    for i in start..end {
        let curr = seq
            .get(Position::new(i).unwrap())
            .unwrap()
            .to_ascii_uppercase();

        if curr == prev {
            hp_len += 1;
        } else {
            if hp_len > 1 {
                res.push(FeatureInterval {
                    start: i - hp_len,
                    stop: i,
                    val: format!("{: >5}{}", hp_len, prev as char),
                });
            }
            hp_len = 1;
            prev = curr;
        }
    }

    if hp_len > 1 {
        res.push(FeatureInterval {
            start: end - hp_len,
            stop: end,
            val: format!("{: >5}{}", hp_len, prev as char),
        });
    }

    res
}

pub fn get_windows(start: usize, end: usize, win_len: usize) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    for i in (start..end).step_by(win_len) {
        res.push(FeatureInterval {
            start: i,
            stop: (i + win_len).min(end),
            val: format!("window_{}", win_len),
        });
    }

    res
}

const BORDER_CONTEXT: usize = 1;

pub fn get_borders(start: usize, end: usize, win_len: usize) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    for i in (start..end).step_by(win_len).skip(1) {
        res.push(FeatureInterval {
            start: (i - 1 - BORDER_CONTEXT).max(start),
            stop: (i + BORDER_CONTEXT + 1).min(end),
            val: format!("border_{}", win_len),
        });
    }

    res
}

pub fn get_matches(
    seq: &fasta::record::Sequence,
    start: usize,
    end: usize,
    seqs: &[String],
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    let mut i = start;
    while i < end {
        for s in seqs {
            // convert to zero-indexed
            if seq.as_ref()[i - 1..(i - 1 + s.len()).min(end - 1)]
                .iter()
                .map(|c| c.to_ascii_uppercase())
                .eq(s.bytes())
            {
                res.push(FeatureInterval {
                    start: i,
                    stop: i + s.len(),
                    val: s.to_owned(),
                });
                i += s.len();
                continue;
            }
        }
        i += 1;
    }

    res
}
