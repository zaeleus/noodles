use super::{Bin, Chunk};

/// A tabix index reference sequence bin builder.
#[derive(Debug, Default)]
pub struct Builder {
    id: u32,
    chunks: Vec<Chunk>,
}

impl Builder {
    /// Set a bin ID.
    pub fn set_id(&mut self, id: u32) -> &mut Self {
        self.id = id;
        self
    }

    /// Adds or merges a chunk.
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

    /// Builds a tabix index reference sequence bin.
    pub fn build(self) -> Bin {
        Bin::new(self.id, self.chunks)
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
