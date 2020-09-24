use noodles_bgzf as bgzf;

/// A tabix index reference sequence bin chunk.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Chunk {
    start: bgzf::VirtualPosition,
    end: bgzf::VirtualPosition,
}

impl Chunk {
    /// Creates a new chunk.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_tabix::index::reference_sequence::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// ```
    pub fn new(start: bgzf::VirtualPosition, end: bgzf::VirtualPosition) -> Self {
        Self { start, end }
    }

    /// Returns the start of the chunk as a virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_tabix::index::reference_sequence::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// assert_eq!(chunk.start(), bgzf::VirtualPosition::from(8));
    /// ```
    pub fn start(&self) -> bgzf::VirtualPosition {
        self.start
    }

    /// Returns the end of the chunk as a virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_tabix::index::reference_sequence::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// assert_eq!(chunk.end(), bgzf::VirtualPosition::from(13));
    /// ```
    pub fn end(&self) -> bgzf::VirtualPosition {
        self.end
    }
}
