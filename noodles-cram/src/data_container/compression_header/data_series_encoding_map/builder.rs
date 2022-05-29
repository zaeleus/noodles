use std::{error, fmt};

use crate::data_container::compression_header::{
    encoding::codec::{Byte, ByteArray, Integer},
    Encoding,
};

use super::DataSeriesEncodingMap;

#[derive(Default)]
pub struct Builder {
    bam_bit_flags_encoding: Option<Encoding<Integer>>,
    cram_bit_flags_encoding: Option<Encoding<Integer>>,
    reference_id_encoding: Option<Encoding<Integer>>,
    read_lengths_encoding: Option<Encoding<Integer>>,
    in_seq_positions_encoding: Option<Encoding<Integer>>,
    read_groups_encoding: Option<Encoding<Integer>>,
    read_names_encoding: Option<Encoding<ByteArray>>,
    next_mate_bit_flags_encoding: Option<Encoding<Integer>>,
    next_fragment_reference_sequence_id_encoding: Option<Encoding<Integer>>,
    next_mate_alignment_start_encoding: Option<Encoding<Integer>>,
    template_size_encoding: Option<Encoding<Integer>>,
    distance_to_next_fragment_encoding: Option<Encoding<Integer>>,
    tag_ids_encoding: Option<Encoding<Integer>>,
    number_of_read_features_encoding: Option<Encoding<Integer>>,
    read_features_codes_encoding: Option<Encoding<Byte>>,
    in_read_positions_encoding: Option<Encoding<Integer>>,
    deletion_lengths_encoding: Option<Encoding<Integer>>,
    stretches_of_bases_encoding: Option<Encoding<ByteArray>>,
    stretches_of_quality_scores_encoding: Option<Encoding<ByteArray>>,
    base_substitution_codes_encoding: Option<Encoding<Byte>>,
    insertion_encoding: Option<Encoding<ByteArray>>,
    reference_skip_length_encoding: Option<Encoding<Integer>>,
    padding_encoding: Option<Encoding<Integer>>,
    hard_clip_encoding: Option<Encoding<Integer>>,
    soft_clip_encoding: Option<Encoding<ByteArray>>,
    mapping_qualities_encoding: Option<Encoding<Integer>>,
    bases_encoding: Option<Encoding<Byte>>,
    quality_scores_encoding: Option<Encoding<Byte>>,
}

impl Builder {
    pub fn set_bam_bit_flags_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.bam_bit_flags_encoding = Some(encoding);
        self
    }

    pub fn set_cram_bit_flags_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.cram_bit_flags_encoding = Some(encoding);
        self
    }

    pub fn set_reference_id_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.reference_id_encoding = Some(encoding);
        self
    }

    pub fn set_read_lengths_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.read_lengths_encoding = Some(encoding);
        self
    }

    pub fn set_in_seq_positions_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.in_seq_positions_encoding = Some(encoding);
        self
    }

    pub fn set_read_groups_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.read_groups_encoding = Some(encoding);
        self
    }

    pub fn set_read_names_encoding(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.read_names_encoding = Some(encoding);
        self
    }

    pub fn set_next_mate_bit_flags_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.next_mate_bit_flags_encoding = Some(encoding);
        self
    }

    pub fn set_next_fragment_reference_sequence_id_encoding(
        mut self,
        encoding: Encoding<Integer>,
    ) -> Self {
        self.next_fragment_reference_sequence_id_encoding = Some(encoding);
        self
    }

    pub fn set_next_mate_alignment_start_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.next_mate_alignment_start_encoding = Some(encoding);
        self
    }

    pub fn set_template_size_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.template_size_encoding = Some(encoding);
        self
    }

    pub fn set_distance_to_next_fragment_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.distance_to_next_fragment_encoding = Some(encoding);
        self
    }

    pub fn set_tag_ids_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.tag_ids_encoding = Some(encoding);
        self
    }

    pub fn set_number_of_read_features_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.number_of_read_features_encoding = Some(encoding);
        self
    }

    pub fn set_read_features_codes_encoding(mut self, encoding: Encoding<Byte>) -> Self {
        self.read_features_codes_encoding = Some(encoding);
        self
    }

    pub fn set_in_read_positions_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.in_read_positions_encoding = Some(encoding);
        self
    }

    pub fn set_deletion_lengths_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.deletion_lengths_encoding = Some(encoding);
        self
    }

    pub fn set_stretches_of_bases_encoding(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.stretches_of_bases_encoding = Some(encoding);
        self
    }

    pub fn set_stretches_of_quality_scores_encoding(
        mut self,
        encoding: Encoding<ByteArray>,
    ) -> Self {
        self.stretches_of_quality_scores_encoding = Some(encoding);
        self
    }

    pub fn set_base_substitution_codes_encoding(mut self, encoding: Encoding<Byte>) -> Self {
        self.base_substitution_codes_encoding = Some(encoding);
        self
    }

    pub fn set_insertion_encoding(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.insertion_encoding = Some(encoding);
        self
    }

    pub fn set_reference_skip_length_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.reference_skip_length_encoding = Some(encoding);
        self
    }

    pub fn set_padding_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.padding_encoding = Some(encoding);
        self
    }

    pub fn set_hard_clip_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.hard_clip_encoding = Some(encoding);
        self
    }

    pub fn set_soft_clip_encoding(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.soft_clip_encoding = Some(encoding);
        self
    }

    pub fn set_mapping_qualities_encoding(mut self, encoding: Encoding<Integer>) -> Self {
        self.mapping_qualities_encoding = Some(encoding);
        self
    }

    pub fn set_bases_encoding(mut self, encoding: Encoding<Byte>) -> Self {
        self.bases_encoding = Some(encoding);
        self
    }

    pub fn set_quality_scores_encoding(mut self, encoding: Encoding<Byte>) -> Self {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.bam_bit_flags_encoding.is_none());
        assert!(builder.cram_bit_flags_encoding.is_none());
        assert!(builder.reference_id_encoding.is_none());
        assert!(builder.read_lengths_encoding.is_none());
        assert!(builder.in_seq_positions_encoding.is_none());
        assert!(builder.read_groups_encoding.is_none());
        assert!(builder.read_names_encoding.is_none());
        assert!(builder.next_mate_bit_flags_encoding.is_none());
        assert!(builder
            .next_fragment_reference_sequence_id_encoding
            .is_none());
        assert!(builder.next_mate_alignment_start_encoding.is_none());
        assert!(builder.template_size_encoding.is_none());
        assert!(builder.distance_to_next_fragment_encoding.is_none());
        assert!(builder.tag_ids_encoding.is_none());
        assert!(builder.number_of_read_features_encoding.is_none());
        assert!(builder.read_features_codes_encoding.is_none());
        assert!(builder.in_read_positions_encoding.is_none());
        assert!(builder.deletion_lengths_encoding.is_none());
        assert!(builder.stretches_of_bases_encoding.is_none());
        assert!(builder.stretches_of_quality_scores_encoding.is_none());
        assert!(builder.base_substitution_codes_encoding.is_none());
        assert!(builder.insertion_encoding.is_none());
        assert!(builder.reference_skip_length_encoding.is_none());
        assert!(builder.padding_encoding.is_none());
        assert!(builder.hard_clip_encoding.is_none());
        assert!(builder.soft_clip_encoding.is_none());
        assert!(builder.mapping_qualities_encoding.is_none());
        assert!(builder.bases_encoding.is_none());
        assert!(builder.quality_scores_encoding.is_none());
    }
}
