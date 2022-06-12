//! CSI reference sequence bin and fields.

mod chunk;

pub use self::chunk::Chunk;

use noodles_bgzf as bgzf;

pub(crate) const METADATA_CHUNK_COUNT: u32 = 2;

/// A CSI reference sequence bin.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bin {
    id: usize,
    loffset: bgzf::VirtualPosition,
    chunks: Vec<Chunk>,
}

impl Bin {
    /// Calculates the maximum bin ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::index::reference_sequence::Bin;
    /// assert_eq!(Bin::max_id(5), 37449);
    /// ```
    pub fn max_id(depth: u8) -> usize {
        bin_limit(depth) as usize
    }

    /// Calculates the metadata bin ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::index::reference_sequence::Bin;
    /// assert_eq!(Bin::metadata_id(5), 37450);
    /// ```
    pub fn metadata_id(depth: u8) -> usize {
        Self::max_id(depth) + 1
    }

    /// Creates a new bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, bgzf::VirtualPosition::default(), Vec::new());
    /// ```
    pub fn new(id: usize, loffset: bgzf::VirtualPosition, chunks: Vec<Chunk>) -> Self {
        Self {
            id,
            loffset,
            chunks,
        }
    }

    /// Returns the bin ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, bgzf::VirtualPosition::default(), Vec::new());
    /// assert_eq!(bin.id(), 10946);
    /// ```
    pub fn id(&self) -> usize {
        self.id
    }

    /// Returns the start virtual position of the first overlapping record in this bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, bgzf::VirtualPosition::default(), Vec::new());
    /// assert_eq!(bin.loffset(), bgzf::VirtualPosition::default());
    /// ```
    pub fn loffset(&self) -> bgzf::VirtualPosition {
        self.loffset
    }

    /// Returns the list of chunks in the bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, bgzf::VirtualPosition::default(), Vec::new());
    /// assert!(bin.chunks().is_empty());
    /// ```
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }
}

// `CSIv1.pdf` (2020-07-21)
fn bin_limit(depth: u8) -> i32 {
    assert!(depth <= 10);
    (1 << ((depth + 1) * 3)) / 7
}
