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
        if let Some(last_chunk) = self.chunks.last_mut()
            && chunk.start() <= last_chunk.end()
        {
            *last_chunk = Chunk::new(last_chunk.start(), chunk.end());
            return;
        }

        self.chunks.push(chunk);
    }

    // Sorts the bin's chunks by their start position.
    pub(crate) fn sort_chunks(&mut self) {
        self.chunks.sort_unstable_by_key(|chunk| chunk.start());
    }

    // Returns the span of the bin's chunks in compressed (BGZF block) offset units, i.e., the
    // distance between the first chunk's start block and the last chunk's end block. Assumes the
    // chunks are sorted by start position.
    pub(crate) fn compressed_span(&self) -> u64 {
        match (self.chunks.first(), self.chunks.last()) {
            (Some(first), Some(last)) => last.end().compressed() - first.start().compressed(),
            _ => 0,
        }
    }

    // Consumes the bin, returning its chunks.
    pub(crate) fn into_chunks(self) -> Vec<Chunk> {
        self.chunks
    }

    // Appends chunks to the bin.
    pub(crate) fn extend_chunks(&mut self, chunks: Vec<Chunk>) {
        self.chunks.extend(chunks);
    }

    // Merges chunks that share a BGZF block into single chunks. Assumes the chunks are sorted by
    // start position. Matches the second pass of htslib's `compress_binning`.
    pub(crate) fn merge_chunks_in_same_block(&mut self) {
        if self.chunks.len() < 2 {
            return;
        }

        let mut m = 0;

        for l in 1..self.chunks.len() {
            if self.chunks[m].end().compressed() >= self.chunks[l].start().compressed() {
                let end = self.chunks[m].end().max(self.chunks[l].end());
                self.chunks[m] = Chunk::new(self.chunks[m].start(), end);
            } else {
                m += 1;
                self.chunks[m] = self.chunks[l];
            }
        }

        self.chunks.truncate(m + 1);
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

    // Builds a virtual position from a compressed (block) offset and an uncompressed offset.
    fn vp(compressed: u64, uncompressed: u16) -> bgzf::VirtualPosition {
        bgzf::VirtualPosition::try_from((compressed, uncompressed)).unwrap()
    }

    #[test]
    fn test_sort_chunks() {
        let mut bin = Bin::new(vec![
            Chunk::new(vp(21, 0), vp(34, 0)),
            Chunk::new(vp(5, 0), vp(13, 0)),
            Chunk::new(vp(8, 0), vp(21, 0)),
        ]);

        bin.sort_chunks();

        assert_eq!(
            bin.chunks,
            [
                Chunk::new(vp(5, 0), vp(13, 0)),
                Chunk::new(vp(8, 0), vp(21, 0)),
                Chunk::new(vp(21, 0), vp(34, 0)),
            ]
        );
    }

    #[test]
    fn test_compressed_span() {
        // last chunk's end block (7) - first chunk's start block (2)
        let bin = Bin::new(vec![
            Chunk::new(vp(2, 0), vp(2, 500)),
            Chunk::new(vp(3, 0), vp(7, 100)),
        ]);
        assert_eq!(bin.compressed_span(), 5);

        // a single chunk contained in one block
        let bin = Bin::new(vec![Chunk::new(vp(4, 10), vp(4, 900))]);
        assert_eq!(bin.compressed_span(), 0);
    }

    #[test]
    fn test_merge_chunks_in_same_block() {
        // two chunks within the same block are merged
        let mut bin = Bin::new(vec![
            Chunk::new(vp(0, 0), vp(0, 50)),
            Chunk::new(vp(0, 50), vp(0, 100)),
        ]);
        bin.merge_chunks_in_same_block();
        assert_eq!(bin.chunks, [Chunk::new(vp(0, 0), vp(0, 100))]);

        // chunks sharing a block across a boundary are merged, taking the max end
        let mut bin = Bin::new(vec![
            Chunk::new(vp(0, 0), vp(1, 10)),
            Chunk::new(vp(1, 20), vp(1, 30)),
        ]);
        bin.merge_chunks_in_same_block();
        assert_eq!(bin.chunks, [Chunk::new(vp(0, 0), vp(1, 30))]);

        // chunks in different blocks are kept separate
        let mut bin = Bin::new(vec![
            Chunk::new(vp(0, 0), vp(0, 50)),
            Chunk::new(vp(1, 0), vp(1, 50)),
        ]);
        bin.merge_chunks_in_same_block();
        assert_eq!(
            bin.chunks,
            [
                Chunk::new(vp(0, 0), vp(0, 50)),
                Chunk::new(vp(1, 0), vp(1, 50))
            ]
        );
    }
}
