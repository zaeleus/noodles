use std::io;

use noodles_core::Position;
use noodles_csi::{
    self as csi,
    binning_index::index::{
        header::ReferenceSequenceNames,
        reference_sequence::{bin::Chunk, index::LinearIndex},
        Header,
    },
};

use crate::Index;

/// A tabix indexer.
#[derive(Debug, Default)]
pub struct Indexer {
    header: Header,
    reference_sequence_names: ReferenceSequenceNames,
    indexer: csi::binning_index::Indexer<LinearIndex>,
}

impl Indexer {
    /// Sets an index header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::binning_index::index::Header;
    /// use noodles_tabix::index::Indexer;
    /// let builder = Indexer::default().set_header(Header::default());
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
    /// use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    /// use noodles_tabix::index::Indexer;
    ///
    /// let mut indexer = Indexer::default();
    ///
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    ///
    /// indexer.add_record("sq0", start, end, Chunk::new(
    ///     bgzf::VirtualPosition::from(144),
    ///     bgzf::VirtualPosition::from(233),
    /// ))?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
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
        self.indexer.add_record(alignment_context, chunk)
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

        *self.header.reference_sequence_names_mut() = self.reference_sequence_names;

        self.indexer
            .set_header(self.header)
            .build(reference_sequence_count)
    }
}
