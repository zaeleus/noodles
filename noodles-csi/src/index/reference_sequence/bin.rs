use noodles_bgzf::{self as bgzf, index::Chunk};

/// A CSI reference sequence bin.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bin {
    id: u32,
    loffset: bgzf::VirtualPosition,
    chunks: Vec<Chunk>,
}

impl Bin {
    /// Creates a new bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, bgzf::VirtualPosition::default(), Vec::new());
    /// ```
    pub fn new(id: u32, loffset: bgzf::VirtualPosition, chunks: Vec<Chunk>) -> Self {
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
    pub fn id(&self) -> u32 {
        self.id
    }

    /// Returns the last offset in the linear index.
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
