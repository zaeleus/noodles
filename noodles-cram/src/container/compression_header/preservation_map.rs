//! CRAM container compression header preservation map.

pub(crate) mod key;
pub(crate) mod substitution_matrix;
pub mod tag_sets;

pub(crate) use {key::Key, substitution_matrix::SubstitutionMatrix, tag_sets::TagSets};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PreservationMap {
    pub(crate) records_have_names: bool,
    pub(crate) alignment_starts_are_deltas: bool,
    pub(crate) external_reference_sequence_is_required: bool,
    pub(crate) substitution_matrix: SubstitutionMatrix,
    pub(crate) tag_sets: TagSets,
}

impl PreservationMap {
    pub fn new(
        records_have_names: bool,
        alignment_starts_are_deltas: bool,
        external_reference_sequence_is_required: bool,
        substitution_matrix: SubstitutionMatrix,
        tag_sets: TagSets,
    ) -> Self {
        Self {
            records_have_names,
            alignment_starts_are_deltas,
            external_reference_sequence_is_required,
            substitution_matrix,
            tag_sets,
        }
    }

    pub fn records_have_names(&self) -> bool {
        self.records_have_names
    }

    pub fn alignment_starts_are_deltas(&self) -> bool {
        self.alignment_starts_are_deltas
    }

    pub fn external_reference_sequence_is_required(&self) -> bool {
        self.external_reference_sequence_is_required
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
            records_have_names: true,
            alignment_starts_are_deltas: true,
            external_reference_sequence_is_required: true,
            substitution_matrix: SubstitutionMatrix::default(),
            tag_sets: TagSets::default(),
        }
    }
}
