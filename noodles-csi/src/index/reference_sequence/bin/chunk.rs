use std::ops::Range;

use noodles_bgzf as bgzf;

/// An index reference sequence bin chunk.
///
/// A chunk is a range of virtual positions representing [start, end).
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
    /// use noodles_csi::index::reference_sequence::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// ```
    pub fn new(start: bgzf::VirtualPosition, end: bgzf::VirtualPosition) -> Self {
        Self { start, end }
    }

    /// Returns the chunk start (inclusive) as a virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// assert_eq!(chunk.start(), bgzf::VirtualPosition::from(8));
    /// ```
    pub fn start(&self) -> bgzf::VirtualPosition {
        self.start
    }

    /// Returns the chunk end (exclusive) as a virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::bin::Chunk;
    /// let chunk = Chunk::new(bgzf::VirtualPosition::from(8), bgzf::VirtualPosition::from(13));
    /// assert_eq!(chunk.end(), bgzf::VirtualPosition::from(13));
    /// ```
    pub fn end(&self) -> bgzf::VirtualPosition {
        self.end
    }
}

impl From<Range<bgzf::VirtualPosition>> for Chunk {
    fn from(range: Range<bgzf::VirtualPosition>) -> Self {
        Self::new(range.start, range.end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_range_bgzf_virtual_position_for_chunk() {
        let start = bgzf::VirtualPosition::from(8);
        let end = bgzf::VirtualPosition::from(13);
        assert_eq!(Chunk::from(start..end), Chunk::new(start, end));
    }
}
