//! CRAM container compression header preservation map.

pub(crate) mod key;
pub(crate) mod substitution_matrix;
pub mod tag_sets;

pub(crate) use {key::Key, substitution_matrix::SubstitutionMatrix, tag_sets::TagSets};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PreservationMap {
    pub(crate) read_names_included: bool,
    pub(crate) ap_data_series_delta: bool,
    pub(crate) is_reference_required: bool,
    pub(crate) substitution_matrix: SubstitutionMatrix,
    pub(crate) tag_sets: TagSets,
}

impl PreservationMap {
    pub fn new(
        read_names_included: bool,
        ap_data_series_delta: bool,
        is_reference_required: bool,
        substitution_matrix: SubstitutionMatrix,
        tag_sets: TagSets,
    ) -> Self {
        Self {
            read_names_included,
            ap_data_series_delta,
            is_reference_required,
            substitution_matrix,
            tag_sets,
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

    pub fn tag_sets(&self) -> &TagSets {
        &self.tag_sets
    }
}

impl Default for PreservationMap {
    fn default() -> Self {
        Self {
            read_names_included: true,
            ap_data_series_delta: true,
            is_reference_required: true,
            substitution_matrix: SubstitutionMatrix::default(),
            tag_sets: TagSets::default(),
        }
    }
}
