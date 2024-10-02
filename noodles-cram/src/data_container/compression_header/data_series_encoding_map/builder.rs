use std::{error, fmt};

use crate::data_container::compression_header::{
    encoding::codec::{Byte, ByteArray, Integer},
    Encoding,
};

use super::DataSeriesEncodingMap;

#[derive(Default)]
pub struct Builder {
    bam_bit_flags: Option<Encoding<Integer>>,
    cram_bit_flags: Option<Encoding<Integer>>,
    reference_id: Option<Encoding<Integer>>,
    read_lengths: Option<Encoding<Integer>>,
    in_seq_positions: Option<Encoding<Integer>>,
    read_groups: Option<Encoding<Integer>>,
    read_names: Option<Encoding<ByteArray>>,
    next_mate_bit_flags: Option<Encoding<Integer>>,
    next_fragment_reference_sequence_id: Option<Encoding<Integer>>,
    next_mate_alignment_start: Option<Encoding<Integer>>,
    template_size: Option<Encoding<Integer>>,
    distance_to_next_fragment: Option<Encoding<Integer>>,
    tag_ids: Option<Encoding<Integer>>,
    number_of_read_features: Option<Encoding<Integer>>,
    read_features_codes: Option<Encoding<Byte>>,
    in_read_positions: Option<Encoding<Integer>>,
    deletion_lengths: Option<Encoding<Integer>>,
    stretches_of_bases: Option<Encoding<ByteArray>>,
    stretches_of_quality_scores: Option<Encoding<ByteArray>>,
    base_substitution_codes: Option<Encoding<Byte>>,
    insertion: Option<Encoding<ByteArray>>,
    reference_skip_length: Option<Encoding<Integer>>,
    padding: Option<Encoding<Integer>>,
    hard_clip: Option<Encoding<Integer>>,
    soft_clip: Option<Encoding<ByteArray>>,
    mapping_qualities: Option<Encoding<Integer>>,
    bases: Option<Encoding<Byte>>,
    quality_scores: Option<Encoding<Byte>>,
}

impl Builder {
    pub fn set_bam_bit_flags(mut self, encoding: Encoding<Integer>) -> Self {
        self.bam_bit_flags = Some(encoding);
        self
    }

    pub fn set_cram_bit_flags(mut self, encoding: Encoding<Integer>) -> Self {
        self.cram_bit_flags = Some(encoding);
        self
    }

    pub fn set_reference_id(mut self, encoding: Encoding<Integer>) -> Self {
        self.reference_id = Some(encoding);
        self
    }

    pub fn set_read_lengths(mut self, encoding: Encoding<Integer>) -> Self {
        self.read_lengths = Some(encoding);
        self
    }

    pub fn set_in_seq_positions(mut self, encoding: Encoding<Integer>) -> Self {
        self.in_seq_positions = Some(encoding);
        self
    }

    pub fn set_read_groups(mut self, encoding: Encoding<Integer>) -> Self {
        self.read_groups = Some(encoding);
        self
    }

    pub fn set_read_names(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.read_names = Some(encoding);
        self
    }

    pub fn set_next_mate_bit_flags(mut self, encoding: Encoding<Integer>) -> Self {
        self.next_mate_bit_flags = Some(encoding);
        self
    }

    pub fn set_next_fragment_reference_sequence_id(mut self, encoding: Encoding<Integer>) -> Self {
        self.next_fragment_reference_sequence_id = Some(encoding);
        self
    }

    pub fn set_next_mate_alignment_start(mut self, encoding: Encoding<Integer>) -> Self {
        self.next_mate_alignment_start = Some(encoding);
        self
    }

    pub fn set_template_size(mut self, encoding: Encoding<Integer>) -> Self {
        self.template_size = Some(encoding);
        self
    }

    pub fn set_distance_to_next_fragment(mut self, encoding: Encoding<Integer>) -> Self {
        self.distance_to_next_fragment = Some(encoding);
        self
    }

    pub fn set_tag_ids(mut self, encoding: Encoding<Integer>) -> Self {
        self.tag_ids = Some(encoding);
        self
    }

    pub fn set_number_of_read_features(mut self, encoding: Encoding<Integer>) -> Self {
        self.number_of_read_features = Some(encoding);
        self
    }

    pub fn set_read_features_codes(mut self, encoding: Encoding<Byte>) -> Self {
        self.read_features_codes = Some(encoding);
        self
    }

    pub fn set_in_read_positions(mut self, encoding: Encoding<Integer>) -> Self {
        self.in_read_positions = Some(encoding);
        self
    }

    pub fn set_deletion_lengths(mut self, encoding: Encoding<Integer>) -> Self {
        self.deletion_lengths = Some(encoding);
        self
    }

    pub fn set_stretches_of_bases(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.stretches_of_bases = Some(encoding);
        self
    }

    pub fn set_stretches_of_quality_scores(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.stretches_of_quality_scores = Some(encoding);
        self
    }

    pub fn set_base_substitution_codes(mut self, encoding: Encoding<Byte>) -> Self {
        self.base_substitution_codes = Some(encoding);
        self
    }

    pub fn set_insertion(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.insertion = Some(encoding);
        self
    }

    pub fn set_reference_skip_length(mut self, encoding: Encoding<Integer>) -> Self {
        self.reference_skip_length = Some(encoding);
        self
    }

    pub fn set_padding(mut self, encoding: Encoding<Integer>) -> Self {
        self.padding = Some(encoding);
        self
    }

    pub fn set_hard_clip(mut self, encoding: Encoding<Integer>) -> Self {
        self.hard_clip = Some(encoding);
        self
    }

    pub fn set_soft_clip(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.soft_clip = Some(encoding);
        self
    }

    pub fn set_mapping_qualities(mut self, encoding: Encoding<Integer>) -> Self {
        self.mapping_qualities = Some(encoding);
        self
    }

    pub fn set_bases(mut self, encoding: Encoding<Byte>) -> Self {
        self.bases = Some(encoding);
        self
    }

    pub fn set_quality_scores(mut self, encoding: Encoding<Byte>) -> Self {
        self.quality_scores = Some(encoding);
        self
    }

    pub(crate) fn build(self) -> Result<DataSeriesEncodingMap, BuildError> {
        Ok(DataSeriesEncodingMap {
            bam_bit_flags: self.bam_bit_flags.ok_or(BuildError::MissingBamBitFlags)?,
            cram_bit_flags: self.cram_bit_flags.ok_or(BuildError::MissingCramBitFlags)?,
            reference_id: self.reference_id,
            read_lengths: self.read_lengths.ok_or(BuildError::MissingReadLengths)?,
            in_seq_positions: self
                .in_seq_positions
                .ok_or(BuildError::MissingInSeqPositions)?,
            read_groups: self.read_groups.ok_or(BuildError::MissingReadGroups)?,
            read_names: self.read_names,
            next_mate_bit_flags: self.next_mate_bit_flags,
            next_fragment_reference_sequence_id: self.next_fragment_reference_sequence_id,
            next_mate_alignment_start: self.next_mate_alignment_start,
            template_size: self.template_size,
            distance_to_next_fragment: self.distance_to_next_fragment,
            tag_ids: self.tag_ids.ok_or(BuildError::MissingTagIds)?,
            number_of_read_features: self.number_of_read_features,
            read_features_codes: self.read_features_codes,
            in_read_positions: self.in_read_positions,
            deletion_lengths: self.deletion_lengths,
            stretches_of_bases: self.stretches_of_bases,
            stretches_of_quality_scores: self.stretches_of_quality_scores,
            base_substitution_codes: self.base_substitution_codes,
            insertion: self.insertion,
            reference_skip_length: self.reference_skip_length,
            padding: self.padding,
            hard_clip: self.hard_clip,
            soft_clip: self.soft_clip,
            mapping_qualities: self.mapping_qualities,
            bases: self.bases,
            quality_scores: self.quality_scores,
        })
    }
}

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    MissingBamBitFlags,
    MissingCramBitFlags,
    MissingReadLengths,
    MissingInSeqPositions,
    MissingReadGroups,
    MissingTagIds,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingBamBitFlags => f.write_str("missing BAM bit flags"),
            Self::MissingCramBitFlags => f.write_str("missing CRAM bit flags"),
            Self::MissingReadLengths => f.write_str("missing read lengths"),
            Self::MissingInSeqPositions => f.write_str("missing in-seq positions"),
            Self::MissingReadGroups => f.write_str("missing read groups"),
            Self::MissingTagIds => f.write_str("missing tag IDs"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.bam_bit_flags.is_none());
        assert!(builder.cram_bit_flags.is_none());
        assert!(builder.reference_id.is_none());
        assert!(builder.read_lengths.is_none());
        assert!(builder.in_seq_positions.is_none());
        assert!(builder.read_groups.is_none());
        assert!(builder.read_names.is_none());
        assert!(builder.next_mate_bit_flags.is_none());
        assert!(builder.next_fragment_reference_sequence_id.is_none());
        assert!(builder.next_mate_alignment_start.is_none());
        assert!(builder.template_size.is_none());
        assert!(builder.distance_to_next_fragment.is_none());
        assert!(builder.tag_ids.is_none());
        assert!(builder.number_of_read_features.is_none());
        assert!(builder.read_features_codes.is_none());
        assert!(builder.in_read_positions.is_none());
        assert!(builder.deletion_lengths.is_none());
        assert!(builder.stretches_of_bases.is_none());
        assert!(builder.stretches_of_quality_scores.is_none());
        assert!(builder.base_substitution_codes.is_none());
        assert!(builder.insertion.is_none());
        assert!(builder.reference_skip_length.is_none());
        assert!(builder.padding.is_none());
        assert!(builder.hard_clip.is_none());
        assert!(builder.soft_clip.is_none());
        assert!(builder.mapping_qualities.is_none());
        assert!(builder.bases.is_none());
        assert!(builder.quality_scores.is_none());
    }
}
