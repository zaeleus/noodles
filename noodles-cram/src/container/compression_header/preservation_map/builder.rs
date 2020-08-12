use super::{
    substitution_matrix, tag_ids_dictionary, PreservationMap, SubstitutionMatrix, TagIdsDictionary,
};

use crate::Record;

pub struct Builder<'a> {
    substitution_matrix_builder: substitution_matrix::Builder<'a>,
    tag_ids_dictionary_builder: tag_ids_dictionary::Builder,
}

impl<'a> Builder<'a> {
    pub fn new(reference_sequence: &'a [u8]) -> Self {
        Self {
            substitution_matrix_builder: SubstitutionMatrix::builder(reference_sequence),
            tag_ids_dictionary_builder: TagIdsDictionary::builder(),
        }
    }

    pub fn update(&mut self, record: &Record) {
        self.substitution_matrix_builder.update(record);
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
