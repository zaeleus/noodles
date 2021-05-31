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
    pub fn new(id: u32, loffset: bgzf::VirtualPosition, chunks: Vec<Chunk>) -> Self {
        Self {
            id,
            loffset,
            chunks,
        }
    }

    /// Returns the bin.
    pub fn id(&self) -> u32 {
        self.id
    }

    /// Returns the list of chunks in the bin.
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }
}
