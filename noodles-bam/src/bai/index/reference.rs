pub mod bin;

pub use self::bin::Bin;

use bit_vec::BitVec;
use noodles_bgzf as bgzf;

const WINDOW_SIZE: u64 = 16384;

#[derive(Debug)]
pub struct Reference {
    bins: Vec<Bin>,
    intervals: Vec<bgzf::VirtualPosition>,
}

impl Reference {
    pub fn new(bins: Vec<Bin>, intervals: Vec<bgzf::VirtualPosition>) -> Self {
        Self { bins, intervals }
    }

    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }

    pub fn intervals(&self) -> &[bgzf::VirtualPosition] {
        &self.intervals
    }

    pub fn query(&self, start: u64, end: u64) -> Vec<&Bin> {
        let region_bins = region_to_bins(start as usize, end as usize);

        let mut query_bins = Vec::new();

        for bin in self.bins() {
            let bin_index = bin.bin() as usize;

            if bin_index < region_bins.len() && region_bins[bin_index] {
                query_bins.push(bin);
            }
        }

        query_bins
    }

    pub fn min_offset(&self, start: u64) -> bgzf::VirtualPosition {
        let i = (start / WINDOW_SIZE) as usize;
        self.intervals.get(i).copied().unwrap_or_default()
    }
}

// § 5.3 C source code for computing bin number and overlapping bins (2020-04-30)
const MAX_BINS: usize = ((1 << 18) - 1) / 7 + 1;

fn region_to_bins(start: usize, end: usize) -> BitVec {
    let ranges = [
        (1 + (start >> 26), 1 + (end >> 26)),
        (9 + (start >> 23), 9 + (end >> 23)),
        (73 + (start >> 20), 73 + (end >> 20)),
        (585 + (start >> 17), 585 + (end >> 17)),
        (4681 + (start >> 14), 4681 + (end >> 14)),
    ];

    let mut bins = BitVec::from_elem(MAX_BINS, false);

    bins.set(0, true);

    for (start, end) in &ranges {
        for k in *start..=*end {
            bins.set(k, true);
        }
    }

    bins
}
