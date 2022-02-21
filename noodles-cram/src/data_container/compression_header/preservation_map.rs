mod builder;
pub mod key;
pub mod substitution_matrix;
mod tag_ids_dictionary;

pub use {
    builder::Builder, key::Key, substitution_matrix::SubstitutionMatrix,
    tag_ids_dictionary::TagIdsDictionary,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PreservationMap {
    read_names_included: bool,
    ap_data_series_delta: bool,
    is_reference_required: bool,
    substitution_matrix: SubstitutionMatrix,
    tag_ids_dictionary: TagIdsDictionary,
}

impl PreservationMap {
    pub fn new(
        read_names_included: bool,
        ap_data_series_delta: bool,
        is_reference_required: bool,
        substitution_matrix: SubstitutionMatrix,
        tag_ids_dictionary: TagIdsDictionary,
    ) -> Self {
        Self {
            read_names_included,
            ap_data_series_delta,
            is_reference_required,
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

    pub fn is_reference_required(&self) -> bool {
        self.is_reference_required
    }

    pub fn substitution_matrix(&self) -> &SubstitutionMatrix {
        &self.substitution_matrix
    }

    pub fn tag_ids_dictionary(&self) -> &TagIdsDictionary {
        &self.tag_ids_dictionary
    }
}
