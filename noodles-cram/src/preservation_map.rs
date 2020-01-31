pub type SubstitutionMatrix = [u8; 5];
pub type TagIdsDictionary = Vec<u8>;

pub struct PreservationMap {
    read_names_included: bool,
    ap_data_series_delta: bool,
    reference_required: bool,
    substitution_matrix: SubstitutionMatrix,
    tag_ids_dictionary: TagIdsDictionary,
}

impl PreservationMap {
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

    pub fn read_names_included_mut(&mut self) -> &mut bool {
        &mut self.read_names_included
    }

    pub fn ap_data_series_delta(&self) -> bool {
        self.ap_data_series_delta
    }

    pub fn ap_data_series_delta_mut(&mut self) -> &mut bool {
        &mut self.ap_data_series_delta
    }

    pub fn reference_required(&self) -> bool {
        self.reference_required
    }

    pub fn reference_required_mut(&mut self) -> &mut bool {
        &mut self.reference_required
    }

    pub fn substitution_matrix(&self) -> &SubstitutionMatrix {
        &self.substitution_matrix
    }

    pub fn substitution_matrix_mut(&mut self) -> &mut SubstitutionMatrix {
        &mut self.substitution_matrix
    }

    pub fn tag_ids_dictionary(&self) -> &TagIdsDictionary {
        &self.tag_ids_dictionary
    }

    pub fn tag_ids_dictionary_mut(&mut self) -> &mut TagIdsDictionary {
        &mut self.tag_ids_dictionary
    }
}

impl Default for PreservationMap {
    fn default() -> Self {
        Self::new(true, true, true, [0, 0, 0, 0, 0], Vec::new())
    }
}
