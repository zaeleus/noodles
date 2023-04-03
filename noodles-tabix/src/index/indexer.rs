use noodles_core::Position;
use noodles_csi::index::{
    header::ReferenceSequenceNames,
    reference_sequence::{self, bin::Chunk},
    Header,
};

use super::Index;

/// A tabix indexer.
#[derive(Debug, Default)]
pub struct Indexer {
    header: Header,
    current_reference_sequence_name: String,
    reference_sequence_names: ReferenceSequenceNames,
    reference_sequence_builders: Vec<reference_sequence::Builder>,
}

impl Indexer {
    /// Sets an index header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use noodles_tabix as tabix;
    /// let builder = tabix::Index::indexer().set_header(csi::index::Header::default());
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
    /// let mut indexer = tabix::Index::indexer();
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
        use super::{DEPTH, MIN_SHIFT};

        if reference_sequence_name != self.current_reference_sequence_name {
            self.reference_sequence_builders
                .push(reference_sequence::Builder::default());

            self.current_reference_sequence_name = reference_sequence_name.into();

            self.reference_sequence_names
                .insert(reference_sequence_name.into());
        }

        let reference_sequence_builder = self
            .reference_sequence_builders
            .last_mut()
            .expect("reference_sequence_builders cannot be empty");

        reference_sequence_builder.add_record(MIN_SHIFT, DEPTH, start, end, true, chunk);
    }

    /// Builds a tabix index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let indexer = tabix::Index::indexer();
    /// let index = indexer.build();
    /// ```
    pub fn build(self) -> Index {
        let reference_sequences = self
            .reference_sequence_builders
            .into_iter()
            .map(|b| b.build())
            .collect();

        Index::builder()
            .set_header(self.header)
            .set_reference_sequence_names(self.reference_sequence_names)
            .set_reference_sequences(reference_sequences)
            .build()
    }
}
