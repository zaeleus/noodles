//! BAM indexer.

pub mod builder;

use std::io;

use noodles_core::Position;
use noodles_csi::{
    self as csi,
    binning_index::index::reference_sequence::{
        bin::Chunk,
        index::{BinnedIndex, LinearIndex},
    },
};

use self::builder::Builder;
use super::Index;

enum Inner {
    Bai(csi::binning_index::Indexer<LinearIndex>),
    Csi(csi::binning_index::Indexer<BinnedIndex>),
}

impl Inner {
    fn add_record(
        &mut self,
        alignment_context: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> io::Result<()> {
        match self {
            Self::Bai(indexer) => indexer.add_record(alignment_context, chunk),
            Self::Csi(indexer) => indexer.add_record(alignment_context, chunk),
        }
    }

    fn build(self, reference_sequence_count: usize) -> Index {
        match self {
            Self::Bai(indexer) => Index::Bai(indexer.build(reference_sequence_count)),
            Self::Csi(indexer) => Index::Csi(indexer.build(reference_sequence_count)),
        }
    }
}

/// A BAM indexer.
pub struct Indexer {
    inner: Inner,
}

impl Indexer {
    /// Creates a BAM indexer builder.
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Adds a record.
    pub fn add_record(
        &mut self,
        alignment_context: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> io::Result<()> {
        self.inner.add_record(alignment_context, chunk)
    }

    /// Builds a BAM index.
    pub fn build(self, reference_sequence_count: usize) -> Index {
        self.inner.build(reference_sequence_count)
    }
}
