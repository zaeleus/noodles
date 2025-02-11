use super::{substitution_matrix, tag_sets, PreservationMap};
use crate::io::writer::{Options, Record};

#[derive(Debug)]
pub struct Builder {
    read_names_included: bool,
    ap_data_series_delta: bool,
    reference_required: bool,
    substitution_matrix_builder: substitution_matrix::Builder,
    tag_ids_dictionary_builder: tag_sets::Builder,
}

impl Builder {
    pub fn apply_options(&mut self, options: &Options) {
        self.read_names_included = options.preserve_read_names;
        self.ap_data_series_delta = options.encode_alignment_start_positions_as_deltas;
    }

    pub fn update(&mut self, record: &Record) {
        self.substitution_matrix_builder.update(record);
        self.tag_ids_dictionary_builder.update(record);
    }

    pub(crate) fn build(self) -> PreservationMap {
        let substitution_matrix = self.substitution_matrix_builder.build();
        let tag_ids_dictionary = self.tag_ids_dictionary_builder.build();

        PreservationMap::new(
            self.read_names_included,
            self.ap_data_series_delta,
            self.reference_required,
            substitution_matrix,
            tag_ids_dictionary,
        )
    }
}

impl Default for Builder {
    // § 8.4 Compression header block (2020-06-22): "The boolean values are optional, defaulting to
    // true when absent, although it is recommended to explicitly set them."
    fn default() -> Self {
        Self {
            read_names_included: true,
            ap_data_series_delta: true,
            reference_required: true,
            substitution_matrix_builder: substitution_matrix::Builder::default(),
            tag_ids_dictionary_builder: tag_sets::Builder::default(),
        }
    }
}
