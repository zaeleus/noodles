//! Coordinate-sorted index (CSI) reference sequence and fields.

mod bin;

pub use self::bin::Bin;

use bit_vec::BitVec;

/// A CSI reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    bins: Vec<Bin>,
}

impl ReferenceSequence {
    /// Creates a CSI reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new());
    /// ```
    pub fn new(bins: Vec<Bin>) -> Self {
        Self { bins }
    }

    /// Returns the list of bins in the reference sequence.
    ///
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new());
    /// assert!(reference_sequence.bins().is_empty());
    /// ```
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }

    pub fn query(&self, min_shift: i32, depth: i32, start: i64, end: i64) -> Vec<&Bin> {
        let max_bin_id = Bin::max_id(depth);
        let mut region_bins = BitVec::from_elem(max_bin_id as usize, false);

        reg2bins(start, end, min_shift, depth, &mut region_bins);

        self.bins()
            .iter()
            .filter(|b| b.id() < max_bin_id && region_bins[b.id() as usize])
            .collect()
    }
}

// `CSIv1.pdf` (2020-07-21)
#[allow(clippy::many_single_char_names)]
fn reg2bins(beg: i64, end: i64, min_shift: i32, depth: i32, bins: &mut BitVec) {
    let mut l = 0;
    let mut t = 0;
    let mut s = min_shift + depth * 3;

    while l <= depth {
        let b = t + (beg >> s);
        let e = t + (end >> s);

        for i in b..=e {
            bins.set(i as usize, true);
        }

        s -= 3;
        t += 1 << (l * 3);
        l += 1;
    }
}
