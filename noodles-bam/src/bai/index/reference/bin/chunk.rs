use noodles_bgzf as bgzf;

/// A chunk in a BAM index bin.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Chunk {
    chunk_beg: bgzf::VirtualPosition,
    chunk_end: bgzf::VirtualPosition,
}

impl Chunk {
    /// Creates a new chunk.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_bam::bai::index::reference::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// ```
    pub fn new(start: bgzf::VirtualPosition, end: bgzf::VirtualPosition) -> Self {
        Self {
            chunk_beg: start,
            chunk_end: end,
        }
    }

    /// The start of the chunk, as a virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_bam::bai::index::reference::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// assert_eq!(chunk.start(), bgzf::VirtualPosition::from(8));
    /// ```
    pub fn start(&self) -> bgzf::VirtualPosition {
        self.chunk_beg
    }

    /// The end of the chunk, as a virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_bam::bai::index::reference::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// assert_eq!(chunk.end(), bgzf::VirtualPosition::from(13));
    /// ```
    pub fn end(&self) -> bgzf::VirtualPosition {
        self.chunk_end
    }
}
