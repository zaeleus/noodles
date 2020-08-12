pub mod key;
pub mod substitution_matrix;
mod tag_ids_dictionary;

pub use {key::Key, substitution_matrix::SubstitutionMatrix, tag_ids_dictionary::TagIdsDictionary};

use crate::Record;

#[derive(Debug)]
pub struct PreservationMap {
    read_names_included: bool,
    ap_data_series_delta: bool,
    reference_required: bool,
    substitution_matrix: SubstitutionMatrix,
    tag_ids_dictionary: TagIdsDictionary,
}

impl PreservationMap {
    pub fn from_records(reference_sequence: &[u8], records: &[Record]) -> Self {
        let mut substitution_matrix_builder = SubstitutionMatrix::builder(reference_sequence);
        let mut tag_ids_dictionary_builder = TagIdsDictionary::builder();

        for record in records {
            substitution_matrix_builder.update(record);
            tag_ids_dictionary_builder.update(record);
        }

        let substitution_matrix = substitution_matrix_builder.build();
        let tag_ids_dictionary = tag_ids_dictionary_builder.build();

        // Read names included, AP data series delta, and reference required all default to `true`.
        // See § 8.4 Compression header block (2020-06-22).
        Self::new(true, true, true, substitution_matrix, tag_ids_dictionary)
    }

    pub fn new(
        read_names_included: bool,
        ap_data_series_delta: bool,
        reference_required: bool,
        substitution_matrix: SubstitutionMatrix,
        tag_ids_dictionary: TagIdsDictionary,
    ) -> Self {
        Self {
            read_names_included,
            ap_data_series_delta,
            reference_required,
            substitution_matrix,
            tag_ids_dictionary,
        }
    }

    pub fn read_names_included(&self) -> bool {
        self.read_names_included
    }

    pub fn ap_data_series_delta(&self) -> bool {
        self.ap_data_series_delta
    }

    pub fn reference_required(&self) -> bool {
        self.reference_required
    }

    pub fn substitution_matrix(&self) -> &SubstitutionMatrix {
        &self.substitution_matrix
    }

    pub fn tag_ids_dictionary(&self) -> &TagIdsDictionary {
        &self.tag_ids_dictionary
    }
}
