use super::{
    reference_sequence::{self, bin::Chunk},
    Index, ReferenceSequence,
};

/// A tabix indexer.
#[derive(Debug, Default)]
pub struct Indexer {
    current_reference_sequence_name: String,
    reference_sequence_names: Vec<String>,
    reference_sequence_builders: Vec<reference_sequence::Builder>,
}

impl Indexer {
    /// Adds a record.
    pub fn add_record(
        &mut self,
        reference_sequence_name: &str,
        start: u32,
        end: u32,
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
    pub fn build(self) -> Index {
        let reference_sequences = self
            .reference_sequence_builders
            .into_iter()
            .map(|b| b.build())
            .collect();

        Index::builder()
            .set_reference_sequence_names(self.reference_sequence_names)
            .set_reference_sequences(reference_sequences)
            .build()
    }
}
