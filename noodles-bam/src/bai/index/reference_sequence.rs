//! BAM index reference sequence and fields.

pub mod bin;

pub use self::bin::Bin;

use bit_vec::BitVec;
use noodles_bgzf as bgzf;

const WINDOW_SIZE: u64 = 16384;

/// A reference sequence in the BAM index.
#[derive(Debug)]
pub struct ReferenceSequence {
    bins: Vec<Bin>,
    intervals: Vec<bgzf::VirtualPosition>,
}

impl ReferenceSequence {
    /// Creates a new BAM index reference seqeuence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// ```
    pub fn new(bins: Vec<Bin>, intervals: Vec<bgzf::VirtualPosition>) -> Self {
        Self { bins, intervals }
    }

    /// Returns the list of bins in this reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// assert!(reference_sequence.bins().is_empty());
    /// ```
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }
    /// Returns a list of 16 kbp intervals that make up the linear index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// assert!(reference_sequence.intervals().is_empty());
    /// ```
    pub fn intervals(&self) -> &[bgzf::VirtualPosition] {
        &self.intervals
    }

    /// Returns a list of bins in this reference seqeunce that intersect the given range.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// let query_bins = reference_sequence.query(8, 13);
    /// assert!(query_bins.is_empty());
    /// ```
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

    /// Finds in minimum start offset in the linear index for a given start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_bam::bai::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// assert_eq!(reference_sequence.min_offset(13), bgzf::VirtualPosition::from(0));
    /// ```
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
