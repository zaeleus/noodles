pub mod bin;

pub use self::bin::Bin;

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

    pub fn min_offset(&self, start: u64) -> bgzf::VirtualPosition {
        let i = (start / WINDOW_SIZE) as usize;
        self.intervals.get(i).copied().unwrap_or_default()
    }
}
