mod builder;
pub mod data_series;

pub use self::{builder::Builder, data_series::DataSeries};

use super::Encoding;

/// A container compression header data series encoding map.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DataSeriesEncodingMap {
    bam_bit_flags_encoding: Encoding,
    cram_bit_flags_encoding: Encoding,
    reference_id_encoding: Option<Encoding>,
    read_lengths_encoding: Encoding,
    in_seq_positions_encoding: Encoding,
    read_groups_encoding: Encoding,
    read_names_encoding: Option<Encoding>,
    next_mate_bit_flags_encoding: Option<Encoding>,
    next_fragment_reference_sequence_id_encoding: Option<Encoding>,
    next_mate_alignment_start_encoding: Option<Encoding>,
    template_size_encoding: Option<Encoding>,
    distance_to_next_fragment_encoding: Option<Encoding>,
    tag_ids_encoding: Encoding,
    number_of_read_features_encoding: Option<Encoding>,
    read_features_codes_encoding: Option<Encoding>,
    in_read_positions_encoding: Option<Encoding>,
    deletion_lengths_encoding: Option<Encoding>,
    stretches_of_bases_encoding: Option<Encoding>,
    stretches_of_quality_scores_encoding: Option<Encoding>,
    base_substitution_codes_encoding: Option<Encoding>,
    insertion_encoding: Option<Encoding>,
    reference_skip_length_encoding: Option<Encoding>,
    padding_encoding: Option<Encoding>,
    hard_clip_encoding: Option<Encoding>,
    soft_clip_encoding: Option<Encoding>,
    mapping_qualities_encoding: Option<Encoding>,
    bases_encoding: Option<Encoding>,
    quality_scores_encoding: Option<Encoding>,
}

impl DataSeriesEncodingMap {
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub fn bam_bit_flags_encoding(&self) -> &Encoding {
        &self.bam_bit_flags_encoding
    }

    pub fn cram_bit_flags_encoding(&self) -> &Encoding {
        &self.cram_bit_flags_encoding
    }

    pub fn reference_id_encoding(&self) -> Option<&Encoding> {
        self.reference_id_encoding.as_ref()
    }

    pub fn read_lengths_encoding(&self) -> &Encoding {
        &self.read_lengths_encoding
    }

    pub fn in_seq_positions_encoding(&self) -> &Encoding {
        &self.in_seq_positions_encoding
    }

    pub fn read_groups_encoding(&self) -> &Encoding {
        &self.read_groups_encoding
    }

    pub fn read_names_encoding(&self) -> Option<&Encoding> {
        self.read_names_encoding.as_ref()
    }

    pub fn next_mate_bit_flags_encoding(&self) -> Option<&Encoding> {
        self.next_mate_bit_flags_encoding.as_ref()
    }

    pub fn next_fragment_reference_sequence_id_encoding(&self) -> Option<&Encoding> {
        self.next_fragment_reference_sequence_id_encoding.as_ref()
    }

    pub fn next_mate_alignment_start_encoding(&self) -> Option<&Encoding> {
        self.next_mate_alignment_start_encoding.as_ref()
    }

    pub fn template_size_encoding(&self) -> Option<&Encoding> {
        self.template_size_encoding.as_ref()
    }

    pub fn distance_to_next_fragment_encoding(&self) -> Option<&Encoding> {
        self.distance_to_next_fragment_encoding.as_ref()
    }

    pub fn tag_ids_encoding(&self) -> &Encoding {
        &self.tag_ids_encoding
    }

    pub fn number_of_read_features_encoding(&self) -> Option<&Encoding> {
        self.number_of_read_features_encoding.as_ref()
    }

    pub fn read_features_codes_encoding(&self) -> Option<&Encoding> {
        self.read_features_codes_encoding.as_ref()
    }

    pub fn in_read_positions_encoding(&self) -> Option<&Encoding> {
        self.in_read_positions_encoding.as_ref()
    }

    pub fn deletion_lengths_encoding(&self) -> Option<&Encoding> {
        self.deletion_lengths_encoding.as_ref()
    }

    pub fn stretches_of_bases_encoding(&self) -> Option<&Encoding> {
        self.stretches_of_bases_encoding.as_ref()
    }

    pub fn stretches_of_quality_scores_encoding(&self) -> Option<&Encoding> {
        self.stretches_of_quality_scores_encoding.as_ref()
    }

    pub fn base_substitution_codes_encoding(&self) -> Option<&Encoding> {
        self.base_substitution_codes_encoding.as_ref()
    }

    pub fn insertion_encoding(&self) -> Option<&Encoding> {
        self.insertion_encoding.as_ref()
    }

    pub fn reference_skip_length_encoding(&self) -> Option<&Encoding> {
        self.reference_skip_length_encoding.as_ref()
    }

    pub fn padding_encoding(&self) -> Option<&Encoding> {
        self.padding_encoding.as_ref()
    }

    pub fn hard_clip_encoding(&self) -> Option<&Encoding> {
        self.hard_clip_encoding.as_ref()
    }

    pub fn soft_clip_encoding(&self) -> Option<&Encoding> {
        self.soft_clip_encoding.as_ref()
    }

    pub fn mapping_qualities_encoding(&self) -> Option<&Encoding> {
        self.mapping_qualities_encoding.as_ref()
    }

    pub fn bases_encoding(&self) -> Option<&Encoding> {
        self.bases_encoding.as_ref()
    }

    pub fn quality_scores_encoding(&self) -> Option<&Encoding> {
        self.quality_scores_encoding.as_ref()
    }
}

impl Default for DataSeriesEncodingMap {
    fn default() -> Self {
        Self {
            bam_bit_flags_encoding: Encoding::External(1),
            cram_bit_flags_encoding: Encoding::External(2),
            reference_id_encoding: Some(Encoding::External(3)),
            read_lengths_encoding: Encoding::External(4),
            in_seq_positions_encoding: Encoding::External(5),
            read_groups_encoding: Encoding::External(6),
            read_names_encoding: Some(Encoding::ByteArrayStop(0x00, 7)),
            next_mate_bit_flags_encoding: Some(Encoding::External(8)),
            next_fragment_reference_sequence_id_encoding: Some(Encoding::External(9)),
            next_mate_alignment_start_encoding: Some(Encoding::External(10)),
            template_size_encoding: Some(Encoding::External(11)),
            distance_to_next_fragment_encoding: Some(Encoding::External(12)),
            tag_ids_encoding: Encoding::External(13),
            number_of_read_features_encoding: Some(Encoding::External(14)),
            read_features_codes_encoding: Some(Encoding::External(15)),
            in_read_positions_encoding: Some(Encoding::External(16)),
            deletion_lengths_encoding: Some(Encoding::External(17)),
            stretches_of_bases_encoding: Some(Encoding::ByteArrayStop(0x00, 18)),
            stretches_of_quality_scores_encoding: Some(Encoding::ByteArrayLen(
                Box::new(Encoding::External(19)),
                Box::new(Encoding::External(19)),
            )),
            base_substitution_codes_encoding: Some(Encoding::External(20)),
            insertion_encoding: Some(Encoding::ByteArrayStop(0x00, 21)),
            reference_skip_length_encoding: Some(Encoding::External(22)),
            padding_encoding: Some(Encoding::External(23)),
            hard_clip_encoding: Some(Encoding::External(24)),
            soft_clip_encoding: Some(Encoding::ByteArrayStop(0x00, 25)),
            mapping_qualities_encoding: Some(Encoding::External(26)),
            bases_encoding: Some(Encoding::External(27)),
            quality_scores_encoding: Some(Encoding::External(28)),
        }
    }
}
