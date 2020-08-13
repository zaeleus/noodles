use super::{substitution_matrix, tag_ids_dictionary, PreservationMap};

use crate::Record;

#[derive(Debug, Default)]
pub struct Builder {
    substitution_matrix_builder: substitution_matrix::Builder,
    tag_ids_dictionary_builder: tag_ids_dictionary::Builder,
}

impl Builder {
    pub fn update(&mut self, reference_sequence: &[u8], record: &Record) {
        self.substitution_matrix_builder
            .update(reference_sequence, record);
        self.tag_ids_dictionary_builder.update(record);
    }

    pub fn build(self) -> PreservationMap {
        let substitution_matrix = self.substitution_matrix_builder.build();
        let tag_ids_dictionary = self.tag_ids_dictionary_builder.build();

        // Read names included, AP data series delta, and reference required all default to `true`.
        // See § 8.4 Compression header block (2020-06-22).
        PreservationMap::new(true, true, true, substitution_matrix, tag_ids_dictionary)
    }
}
