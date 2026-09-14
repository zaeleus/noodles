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
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::index::Indexer;
    /// let builder = Indexer::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Adds a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::index::Indexer;
    /// use noodles_bgzf as bgzf;
    /// use noodles_core::Position;
    /// use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    ///
    /// let mut indexer = Indexer::builder().build()?;
    ///
    /// let reference_sequence_id = 0;
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    /// let is_mapped = true;
    /// let alignment_context = Some((reference_sequence_id, start, end, is_mapped));
    ///
    /// let chunk = Chunk::new(
    ///     bgzf::VirtualPosition::from(144),
    ///     bgzf::VirtualPosition::from(233),
    /// );
    ///
    /// indexer.add_record(alignment_context, chunk)?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn add_record(
        &mut self,
        alignment_context: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> io::Result<()> {
        self.inner.add_record(alignment_context, chunk)
    }

    /// Builds a BAM index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::index::Indexer;
    /// let indexer = Indexer::builder().build()?;
    /// let index = indexer.build(0);
    /// # Ok::<_, noodles_bam::index::indexer::builder::BuildError>(())
    /// ```
    pub fn build(self, reference_sequence_count: usize) -> Index {
        self.inner.build(reference_sequence_count)
    }
}
