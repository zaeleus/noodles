use std::{error, fmt};

use crate::container::compression_header::Encoding;

use super::DataSeriesEncodingMap;

pub struct Builder {
    bam_bit_flags_encoding: Option<Encoding>,
    cram_bit_flags_encoding: Option<Encoding>,
    reference_id_encoding: Option<Encoding>,
    read_lengths_encoding: Option<Encoding>,
    in_seq_positions_encoding: Option<Encoding>,
    read_groups_encoding: Option<Encoding>,
    read_names_encoding: Option<Encoding>,
    next_mate_bit_flags_encoding: Option<Encoding>,
    next_fragment_reference_sequence_id_encoding: Option<Encoding>,
    next_mate_alignment_start_encoding: Option<Encoding>,
    template_size_encoding: Option<Encoding>,
    distance_to_next_fragment_encoding: Option<Encoding>,
    tag_ids_encoding: Option<Encoding>,
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

impl Builder {
    pub fn set_bam_bit_flags_encoding(mut self, encoding: Encoding) -> Self {
        self.bam_bit_flags_encoding = Some(encoding);
        self
    }

    pub fn set_cram_bit_flags_encoding(mut self, encoding: Encoding) -> Self {
        self.cram_bit_flags_encoding = Some(encoding);
        self
    }

    pub fn set_reference_id_encoding(mut self, encoding: Encoding) -> Self {
        self.reference_id_encoding = Some(encoding);
        self
    }

    pub fn set_read_lengths_encoding(mut self, encoding: Encoding) -> Self {
        self.read_lengths_encoding = Some(encoding);
        self
    }

    pub fn set_in_seq_positions_encoding(mut self, encoding: Encoding) -> Self {
        self.in_seq_positions_encoding = Some(encoding);
        self
    }

    pub fn set_read_groups_encoding(mut self, encoding: Encoding) -> Self {
        self.read_groups_encoding = Some(encoding);
        self
    }

    pub fn set_read_names_encoding(mut self, encoding: Encoding) -> Self {
        self.read_names_encoding = Some(encoding);
        self
    }

    pub fn set_next_mate_bit_flags_encoding(mut self, encoding: Encoding) -> Self {
        self.next_mate_bit_flags_encoding = Some(encoding);
        self
    }

    pub fn set_next_fragment_reference_sequence_id_encoding(mut self, encoding: Encoding) -> Self {
        self.next_fragment_reference_sequence_id_encoding = Some(encoding);
        self
    }

    pub fn set_next_mate_alignment_start_encoding(mut self, encoding: Encoding) -> Self {
        self.next_mate_alignment_start_encoding = Some(encoding);
        self
    }

    pub fn set_template_size_encoding(mut self, encoding: Encoding) -> Self {
        self.template_size_encoding = Some(encoding);
        self
    }

    pub fn set_distance_to_next_fragment_encoding(mut self, encoding: Encoding) -> Self {
        self.distance_to_next_fragment_encoding = Some(encoding);
        self
    }

    pub fn set_tag_ids_encoding(mut self, encoding: Encoding) -> Self {
        self.tag_ids_encoding = Some(encoding);
        self
    }

    pub fn set_number_of_read_features_encoding(mut self, encoding: Encoding) -> Self {
        self.number_of_read_features_encoding = Some(encoding);
        self
    }

    pub fn set_read_features_codes_encoding(mut self, encoding: Encoding) -> Self {
        self.read_features_codes_encoding = Some(encoding);
        self
    }

    pub fn set_in_read_positions_encoding(mut self, encoding: Encoding) -> Self {
        self.in_read_positions_encoding = Some(encoding);
        self
    }

    pub fn set_deletion_lengths_encoding(mut self, encoding: Encoding) -> Self {
        self.deletion_lengths_encoding = Some(encoding);
        self
    }

    pub fn set_stretches_of_bases_encoding(mut self, encoding: Encoding) -> Self {
        self.stretches_of_bases_encoding = Some(encoding);
        self
    }

    pub fn set_stretches_of_quality_scores_encoding(mut self, encoding: Encoding) -> Self {
        self.stretches_of_quality_scores_encoding = Some(encoding);
        self
    }

    pub fn set_base_substitution_codes_encoding(mut self, encoding: Encoding) -> Self {
        self.base_substitution_codes_encoding = Some(encoding);
        self
    }

    pub fn set_insertion_encoding(mut self, encoding: Encoding) -> Self {
        self.insertion_encoding = Some(encoding);
        self
    }

    pub fn set_reference_skip_length_encoding(mut self, encoding: Encoding) -> Self {
        self.reference_skip_length_encoding = Some(encoding);
        self
    }

    pub fn set_padding_encoding(mut self, encoding: Encoding) -> Self {
        self.padding_encoding = Some(encoding);
        self
    }

    pub fn set_hard_clip_encoding(mut self, encoding: Encoding) -> Self {
        self.hard_clip_encoding = Some(encoding);
        self
    }

    pub fn set_soft_clip_encoding(mut self, encoding: Encoding) -> Self {
        self.soft_clip_encoding = Some(encoding);
        self
    }

    pub fn set_mapping_qualities_encoding(mut self, encoding: Encoding) -> Self {
        self.mapping_qualities_encoding = Some(encoding);
        self
    }

    pub fn set_bases_encoding(mut self, encoding: Encoding) -> Self {
        self.bases_encoding = Some(encoding);
        self
    }

    pub fn set_quality_scores_encoding(mut self, encoding: Encoding) -> Self {
        self.quality_scores_encoding = Some(encoding);
        self
    }

    pub fn build(self) -> Result<DataSeriesEncodingMap, BuildError> {
        Ok(DataSeriesEncodingMap {
            bam_bit_flags_encoding: self
                .bam_bit_flags_encoding
                .ok_or(BuildError::MissingBamBitFlagsEncoding)?,
            cram_bit_flags_encoding: self
                .cram_bit_flags_encoding
                .ok_or(BuildError::MissingCramBitFlagsEncoding)?,
            reference_id_encoding: self.reference_id_encoding,
            read_lengths_encoding: self
                .read_lengths_encoding
                .ok_or(BuildError::MissingReadLengthsEncoding)?,
            in_seq_positions_encoding: self
                .in_seq_positions_encoding
                .ok_or(BuildError::MissingInSeqPositionsEncoding)?,
            read_groups_encoding: self
                .read_groups_encoding
                .ok_or(BuildError::MissingReadGroupsEncoding)?,
            read_names_encoding: self.read_names_encoding,
            next_mate_bit_flags_encoding: self.next_mate_bit_flags_encoding,
            next_fragment_reference_sequence_id_encoding: self
                .next_fragment_reference_sequence_id_encoding,
            next_mate_alignment_start_encoding: self.next_mate_alignment_start_encoding,
            template_size_encoding: self.template_size_encoding,
            distance_to_next_fragment_encoding: self.distance_to_next_fragment_encoding,
            tag_ids_encoding: self
                .tag_ids_encoding
                .ok_or(BuildError::MissingTagIdsEncoding)?,
            number_of_read_features_encoding: self.number_of_read_features_encoding,
            read_features_codes_encoding: self.read_features_codes_encoding,
            in_read_positions_encoding: self.in_read_positions_encoding,
            deletion_lengths_encoding: self.deletion_lengths_encoding,
            stretches_of_bases_encoding: self.stretches_of_bases_encoding,
            stretches_of_quality_scores_encoding: self.stretches_of_quality_scores_encoding,
            base_substitution_codes_encoding: self.base_substitution_codes_encoding,
            insertion_encoding: self.insertion_encoding,
            reference_skip_length_encoding: self.reference_skip_length_encoding,
            padding_encoding: self.padding_encoding,
            hard_clip_encoding: self.hard_clip_encoding,
            soft_clip_encoding: self.soft_clip_encoding,
            mapping_qualities_encoding: self.mapping_qualities_encoding,
            bases_encoding: self.bases_encoding,
            quality_scores_encoding: self.quality_scores_encoding,
        })
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            bam_bit_flags_encoding: Some(Encoding::External(1)),
            cram_bit_flags_encoding: Some(Encoding::External(2)),
            reference_id_encoding: Some(Encoding::External(3)),
            read_lengths_encoding: Some(Encoding::External(4)),
            in_seq_positions_encoding: Some(Encoding::External(5)),
            read_groups_encoding: Some(Encoding::External(6)),
            read_names_encoding: Some(Encoding::ByteArrayStop(0x00, 7)),
            next_mate_bit_flags_encoding: Some(Encoding::External(8)),
            next_fragment_reference_sequence_id_encoding: Some(Encoding::External(9)),
            next_mate_alignment_start_encoding: Some(Encoding::External(10)),
            template_size_encoding: Some(Encoding::External(11)),
            distance_to_next_fragment_encoding: Some(Encoding::External(12)),
            tag_ids_encoding: Some(Encoding::External(13)),
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

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    MissingBamBitFlagsEncoding,
    MissingCramBitFlagsEncoding,
    MissingReadLengthsEncoding,
    MissingInSeqPositionsEncoding,
    MissingReadGroupsEncoding,
    MissingTagIdsEncoding,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingBamBitFlagsEncoding => f.write_str("missing BAM bit flags encoding"),
            Self::MissingCramBitFlagsEncoding => f.write_str("missing CRAM bit flags encoding"),
            Self::MissingReadLengthsEncoding => f.write_str("missing read lengths encoding"),
            Self::MissingInSeqPositionsEncoding => f.write_str("missing in-seq positions encoding"),
            Self::MissingReadGroupsEncoding => f.write_str("missing read groups encoding"),
            Self::MissingTagIdsEncoding => f.write_str("missing tag IDs encoding"),
        }
    }
}
