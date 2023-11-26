//! Binning index reference sequence bin.

mod chunk;

pub use self::chunk::Chunk;

pub(crate) const METADATA_CHUNK_COUNT: u32 = 2;

/// A binning index reference sequence bin.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bin {
    chunks: Vec<Chunk>,
}

impl Bin {
    /// Calculates the maximum bin ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::index::reference_sequence::Bin;
    /// assert_eq!(Bin::max_id(5), 37449);
    /// ```
    pub const fn max_id(depth: u8) -> usize {
        bin_limit(depth) as usize
    }

    /// Calculates the metadata bin ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::index::reference_sequence::Bin;
    /// assert_eq!(Bin::metadata_id(5), 37450);
    /// ```
    pub const fn metadata_id(depth: u8) -> usize {
        Self::max_id(depth) + 1
    }

    /// Creates a binning index reference sequence bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::index::reference_sequence::Bin;
    /// let bin = Bin::new(Vec::new());
    /// ```
    pub fn new(chunks: Vec<Chunk>) -> Self {
        Self { chunks }
    }

    /// Returns the list of chunks in the bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::binning_index::index::reference_sequence::Bin;
    /// let bin = Bin::new(Vec::new());
    /// assert!(bin.chunks().is_empty());
    /// ```
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }

    /// Adds or merges a chunk.
    pub fn add_chunk(&mut self, chunk: Chunk) {
        if let Some(last_chunk) = self.chunks.last_mut() {
            if chunk.start() <= last_chunk.end() {
                *last_chunk = Chunk::new(last_chunk.start(), chunk.end());
                return;
            }
        }

        self.chunks.push(chunk);
    }
}

// `CSIv1.pdf` (2020-07-21)
const fn bin_limit(depth: u8) -> i32 {
    assert!(depth <= 10);
    (1 << ((depth + 1) * 3)) / 7
}

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;

    use super::*;

    #[test]
    fn test_add_chunk() {
        let mut bin = Bin::new(Vec::new());

        bin.add_chunk(Chunk::new(
            bgzf::VirtualPosition::from(5),
            bgzf::VirtualPosition::from(13),
        ));

        assert_eq!(
            bin.chunks,
            [Chunk::new(
                bgzf::VirtualPosition::from(5),
                bgzf::VirtualPosition::from(13)
            )]
        );

        bin.add_chunk(Chunk::new(
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(21),
        ));

        assert_eq!(
            bin.chunks,
            [Chunk::new(
                bgzf::VirtualPosition::from(5),
                bgzf::VirtualPosition::from(21)
            )]
        );

        bin.add_chunk(Chunk::new(
            bgzf::VirtualPosition::from(34),
            bgzf::VirtualPosition::from(55),
        ));

        assert_eq!(
            bin.chunks,
            [
                Chunk::new(
                    bgzf::VirtualPosition::from(5),
                    bgzf::VirtualPosition::from(21)
                ),
                Chunk::new(
                    bgzf::VirtualPosition::from(34),
                    bgzf::VirtualPosition::from(55)
                )
            ]
        );
    }
}
