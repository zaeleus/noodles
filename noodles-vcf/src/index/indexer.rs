mod builder;

use std::io;

use noodles_core::Position;
use noodles_csi::{
    self as csi,
    binning_index::index::{
        header::ReferenceSequenceNames,
        reference_sequence::{
            bin::Chunk,
            index::{BinnedIndex, LinearIndex},
        },
    },
};

use self::builder::Builder;
use super::Index;

enum Inner {
    Csi(csi::binning_index::Indexer<BinnedIndex>),
    Tabix(csi::binning_index::Indexer<LinearIndex>),
}

impl Inner {
    fn set_header(self, header: csi::binning_index::index::Header) -> Self {
        match self {
            Inner::Csi(indexer) => Inner::Csi(indexer.set_header(header)),
            Inner::Tabix(indexer) => Inner::Tabix(indexer.set_header(header)),
        }
    }

    fn add_record(
        &mut self,
        alignment_context: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> io::Result<()> {
        match self {
            Self::Csi(indexer) => indexer.add_record(alignment_context, chunk),
            Self::Tabix(indexer) => indexer.add_record(alignment_context, chunk),
        }
    }

    fn build(self, reference_sequence_count: usize) -> Index {
        match self {
            Self::Csi(indexer) => Index::Csi(indexer.build(reference_sequence_count)),
            Self::Tabix(indexer) => Index::Tabix(indexer.build(reference_sequence_count)),
        }
    }
}

/// A VCF indexer.
pub struct Indexer {
    header: csi::binning_index::index::Header,
    reference_sequence_names: ReferenceSequenceNames,
    inner: Inner,
}

impl Indexer {
    /// Creates a VCF indexer builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::index::Indexer;
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
    /// use noodles_bgzf as bgzf;
    /// use noodles_core::Position;
    /// use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    /// use noodles_vcf::index::Indexer;
    ///
    /// let mut indexer = Indexer::builder().build();
    ///
    /// let reference_sequence_name = "sq0";
    /// let start = const { Position::new(8).unwrap() };
    /// let end = const { Position::new(13).unwrap() };
    /// let chunk = Chunk::new(
    ///     bgzf::VirtualPosition::from(144),
    ///     bgzf::VirtualPosition::from(233),
    /// );
    ///
    /// indexer.add_record(reference_sequence_name, start, end, chunk)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn add_record(
        &mut self,
        reference_sequence_name: &str,
        start: Position,
        end: Position,
        chunk: Chunk,
    ) -> io::Result<()> {
        let (reference_sequence_id, _) = self
            .reference_sequence_names
            .insert_full(reference_sequence_name.into());

        let alignment_context = Some((reference_sequence_id, start, end, true));

        self.inner.add_record(alignment_context, chunk)
    }

    /// Builds a VCF index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::index::Indexer;
    /// let indexer = Indexer::builder().build();
    /// let index = indexer.build();
    /// ```
    pub fn build(mut self) -> Index {
        let reference_sequence_count = self.reference_sequence_names.len();

        *self.header.reference_sequence_names_mut() = self.reference_sequence_names;

        self.inner
            .set_header(self.header)
            .build(reference_sequence_count)
    }
}
