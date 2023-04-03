use noodles_core::Position;
use noodles_csi::{
    self as csi,
    index::{header::ReferenceSequenceNames, reference_sequence::bin::Chunk, Header},
    Index,
};

/// A tabix indexer.
#[derive(Debug, Default)]
pub struct Indexer {
    header: Header,
    reference_sequence_names: ReferenceSequenceNames,
    indexer: csi::index::Indexer,
}

impl Indexer {
    /// Sets an index header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use noodles_tabix as tabix;
    ///
    /// let builder = tabix::index::Indexer::default()
    ///     .set_header(csi::index::Header::default());
    /// ```
    pub fn set_header(&mut self, header: Header) {
        self.header = header;
    }

    /// Adds a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_core::Position;
    /// use noodles_csi::index::reference_sequence::bin::Chunk;
    /// use noodles_tabix as tabix;
    ///
    /// let mut indexer = tabix::index::Indexer::default();
    ///
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    ///
    /// indexer.add_record("sq0", start, end, Chunk::new(
    ///     bgzf::VirtualPosition::from(144),
    ///     bgzf::VirtualPosition::from(233),
    /// ));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn add_record(
        &mut self,
        reference_sequence_name: &str,
        start: Position,
        end: Position,
        chunk: Chunk,
    ) {
        // TODO: Validate contiguous regions.
        let (reference_sequence_id, _) = self
            .reference_sequence_names
            .insert_full(reference_sequence_name.into());

        let alignment_context = Some((reference_sequence_id, start, end, true));
        self.indexer.add_record(alignment_context, chunk).unwrap();
    }

    /// Builds a tabix index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let index = tabix::index::Indexer::default().build();
    /// ```
    pub fn build(mut self) -> Index {
        let reference_sequence_count = self.reference_sequence_names.len();

        self.header.reference_sequence_names = self.reference_sequence_names;

        self.indexer
            .set_header(self.header)
            .build(reference_sequence_count)
    }
}
