use noodles_csi::index::reference_sequence::bin::Chunk;

use super::{reference_sequence, Header, Index, ReferenceSequence};

/// A tabix indexer.
#[derive(Debug, Default)]
pub struct Indexer {
    header: Header,
    current_reference_sequence_name: String,
    reference_sequence_names: Vec<String>,
    reference_sequence_builders: Vec<reference_sequence::Builder>,
}

impl Indexer {
    /// Sets an index header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let builder = tabix::Index::indexer().set_header(tabix::index::Header::default());
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
    /// use noodles_csi::index::reference_sequence::bin::Chunk;
    /// use noodles_tabix as tabix;
    ///
    /// let mut indexer = tabix::Index::indexer();
    ///
    /// indexer.add_record("sq0", 8, 13, Chunk::new(
    ///     bgzf::VirtualPosition::from(144),
    ///     bgzf::VirtualPosition::from(233),
    /// ));
    /// ```
    pub fn add_record(
        &mut self,
        reference_sequence_name: &str,
        start: i32,
        end: i32,
        chunk: Chunk,
    ) {
        if reference_sequence_name != self.current_reference_sequence_name {
            self.reference_sequence_builders
                .push(ReferenceSequence::builder());

            self.current_reference_sequence_name = reference_sequence_name.into();

            self.reference_sequence_names
                .push(reference_sequence_name.into());
        }

        let reference_sequence_builder = self
            .reference_sequence_builders
            .last_mut()
            .expect("reference_sequence_builders cannot be empty");

        reference_sequence_builder.add_record(start, end, chunk);
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
