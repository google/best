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

use noodles::core::Position;
use noodles::fasta;

use crate::bed::FeatureInterval;

static COMPLEMENT: [u8; 128] = {
    let mut c = [0u8; 128];
    c[b'A' as usize] = b'T';
    c[b'T' as usize] = b'A';
    c[b'C' as usize] = b'G';
    c[b'G' as usize] = b'C';
    c[b'N' as usize] = b'N';
    c
};

/// Find homopolymers in a sequence to use as intervals.
pub fn find_homopolymers(
    seq: &fasta::record::Sequence,
    start: usize,
    end: usize,
    strand_rev: bool,
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
                let c = if strand_rev {
                    COMPLEMENT[prev as usize]
                } else {
                    prev
                };
                res.push(FeatureInterval {
                    start: i - hp_len,
                    stop: i,
                    val: format!("{: >5}{}", hp_len, c as char),
                });
            }
            hp_len = 1;
            prev = curr;
        }
    }

    if hp_len > 1 {
        let c = if strand_rev {
            COMPLEMENT[prev as usize]
        } else {
            prev
        };
        res.push(FeatureInterval {
            start: end - hp_len,
            stop: end,
            val: format!("{: >5}{}", hp_len, c as char),
        });
    }

    res
}

/// Get fixed-length windows as intervals.
pub fn get_windows(
    start: usize,
    end: usize,
    win_len: usize,
    pos: bool,
    strand_rev: bool,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    for i in (start..end).step_by(win_len) {
        let lo;
        let hi;
        if strand_rev {
            hi = end - (i - start);
            lo = hi.saturating_sub(win_len).max(start);
        } else {
            lo = i;
            hi = (i + win_len).min(end);
        };

        if pos {
            res.push(FeatureInterval {
                start: lo,
                stop: hi,
                val: format!("window_{}_pos_{}", win_len, i - start),
            });
        } else {
            res.push(FeatureInterval {
                start: lo,
                stop: hi,
                val: format!("window_{}", win_len),
            });
        }
    }

    res
}

const BORDER_CONTEXT: usize = 1;

/// Get small intervals that represent the region near fixed-width window borders.
pub fn get_borders(
    start: usize,
    end: usize,
    win_len: usize,
    strand_rev: bool,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    for i in (start..end).step_by(win_len).skip(1) {
        let idx;
        if strand_rev {
            idx = end - (i - start);
        } else {
            idx = i;
        };

        res.push(FeatureInterval {
            start: (idx - 1 - BORDER_CONTEXT).max(start),
            stop: (idx + BORDER_CONTEXT + 1).min(end),
            val: format!("border_{}", win_len),
        });
    }

    res
}

/// Get regions that match a sequence as intervals.
pub fn get_matches(
    seq: &fasta::record::Sequence,
    start: usize,
    end: usize,
    s: &str,
    strand_rev: bool,
) -> Vec<FeatureInterval> {
    let mut res = Vec::new();

    for i in start..end {
        // convert to zero-indexed
        let seq_iter = seq.as_ref()[i - 1..(i - 1 + s.len()).min(end - 1)]
            .iter()
            .map(|c| c.to_ascii_uppercase());
        let is_match = if strand_rev {
            seq_iter.eq(s.bytes().rev().map(|c| COMPLEMENT[c as usize]))
        } else {
            seq_iter.eq(s.bytes())
        };
        if is_match {
            res.push(FeatureInterval {
                start: i,
                stop: i + s.len(),
                val: s.to_owned(),
            });
        }
    }

    res
}
