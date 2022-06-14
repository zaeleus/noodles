use noodles_csi::index::reference_sequence::bin::Chunk;

use super::Bin;

/// A BAM index reference sequence bin builder.
#[derive(Debug, Default)]
pub struct Builder {
    id: usize,
    chunks: Vec<Chunk>,
}

impl Builder {
    /// Sets a bin ID.
    pub fn set_id(&mut self, id: usize) -> &mut Self {
        self.id = id;
        self
    }

    /// Adds or merges a chunk.
    ///
    /// If the given chunk overlaps the last chunk, it is merged into the last chunk. For example,
    /// adding [2, 5] and then [3, 8] will produce the list [[2, 8]]. Subsequently adding [13, 21],
    /// the final list will be [[2, 8], [13, 21]].
    ///
    /// See § 5.1.2 Reducing small chunks (2020-07-19).
    pub fn add_chunk(&mut self, chunk: Chunk) -> &mut Self {
        if let Some(last_chunk) = self.chunks.last_mut() {
            if chunk.start() <= last_chunk.end() {
                *last_chunk = Chunk::new(last_chunk.start(), chunk.end());
                return self;
            }
        }

        self.chunks.push(chunk);

        self
    }

    /// Builds a BAM index reference sequence bin.
    pub fn build(self) -> Bin {
        Bin {
            id: self.id,
            chunks: self.chunks,
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;

    use super::*;

    #[test]
    fn test_set_id() {
        let mut builder = Builder::default();
        builder.set_id(13);
        assert_eq!(builder.id, 13);
    }

    #[test]
    fn test_add_chunk() {
        let mut builder = Builder::default();

        assert!(builder.chunks.is_empty());

        builder.add_chunk(Chunk::new(
            bgzf::VirtualPosition::from(5),
            bgzf::VirtualPosition::from(13),
        ));

        assert_eq!(
            builder.chunks,
            [Chunk::new(
                bgzf::VirtualPosition::from(5),
                bgzf::VirtualPosition::from(13)
            )]
        );

        builder.add_chunk(Chunk::new(
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(21),
        ));

        assert_eq!(
            builder.chunks,
            [Chunk::new(
                bgzf::VirtualPosition::from(5),
                bgzf::VirtualPosition::from(21)
            )]
        );

        builder.add_chunk(Chunk::new(
            bgzf::VirtualPosition::from(34),
            bgzf::VirtualPosition::from(55),
        ));

        assert_eq!(
            builder.chunks,
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

    #[test]
    fn test_build() {
        let mut builder = Builder::default();
        builder.set_id(13);

        builder.add_chunk(Chunk::new(
            bgzf::VirtualPosition::from(5),
            bgzf::VirtualPosition::from(13),
        ));

        let bin = builder.build();

        assert_eq!(bin.id(), 13);
        assert_eq!(
            bin.chunks(),
            [Chunk::new(
                bgzf::VirtualPosition::from(5),
                bgzf::VirtualPosition::from(13),
            )]
        )
    }
}
