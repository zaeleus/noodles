use std::{error, fmt};

use crate::container::compression_header::{
    encoding::codec::{Byte, ByteArray, Integer},
    Encoding,
};

use super::DataSeriesEncodings;

#[derive(Default)]
pub struct Builder {
    bam_flags: Option<Encoding<Integer>>,
    cram_flags: Option<Encoding<Integer>>,
    reference_sequence_ids: Option<Encoding<Integer>>,
    read_lengths: Option<Encoding<Integer>>,
    alignment_starts: Option<Encoding<Integer>>,
    read_group_ids: Option<Encoding<Integer>>,
    names: Option<Encoding<ByteArray>>,
    mate_flags: Option<Encoding<Integer>>,
    mate_reference_sequence_ids: Option<Encoding<Integer>>,
    mate_alignment_starts: Option<Encoding<Integer>>,
    template_lengths: Option<Encoding<Integer>>,
    mate_distances: Option<Encoding<Integer>>,
    tag_set_ids: Option<Encoding<Integer>>,
    feature_counts: Option<Encoding<Integer>>,
    feature_codes: Option<Encoding<Byte>>,
    feature_position_deltas: Option<Encoding<Integer>>,
    deletion_lengths: Option<Encoding<Integer>>,
    stretches_of_bases: Option<Encoding<ByteArray>>,
    stretches_of_quality_scores: Option<Encoding<ByteArray>>,
    base_substitution_codes: Option<Encoding<Byte>>,
    insertion_bases: Option<Encoding<ByteArray>>,
    reference_skip_lengths: Option<Encoding<Integer>>,
    padding_lengths: Option<Encoding<Integer>>,
    hard_clip_lengths: Option<Encoding<Integer>>,
    soft_clip_bases: Option<Encoding<ByteArray>>,
    mapping_qualities: Option<Encoding<Integer>>,
    bases: Option<Encoding<Byte>>,
    quality_scores: Option<Encoding<Byte>>,
}

impl Builder {
    pub fn set_bam_flags(mut self, encoding: Encoding<Integer>) -> Self {
        self.bam_flags = Some(encoding);
        self
    }

    pub fn set_cram_flags(mut self, encoding: Encoding<Integer>) -> Self {
        self.cram_flags = Some(encoding);
        self
    }

    pub fn set_reference_sequence_ids(mut self, encoding: Encoding<Integer>) -> Self {
        self.reference_sequence_ids = Some(encoding);
        self
    }

    pub fn set_read_lengths(mut self, encoding: Encoding<Integer>) -> Self {
        self.read_lengths = Some(encoding);
        self
    }

    pub fn set_alignment_starts(mut self, encoding: Encoding<Integer>) -> Self {
        self.alignment_starts = Some(encoding);
        self
    }

    pub fn set_read_group_ids(mut self, encoding: Encoding<Integer>) -> Self {
        self.read_group_ids = Some(encoding);
        self
    }

    pub fn set_names(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.names = Some(encoding);
        self
    }

    pub fn set_mate_flags(mut self, encoding: Encoding<Integer>) -> Self {
        self.mate_flags = Some(encoding);
        self
    }

    pub fn set_mate_reference_sequence_ids(mut self, encoding: Encoding<Integer>) -> Self {
        self.mate_reference_sequence_ids = Some(encoding);
        self
    }

    pub fn set_mate_alignment_starts(mut self, encoding: Encoding<Integer>) -> Self {
        self.mate_alignment_starts = Some(encoding);
        self
    }

    pub fn set_template_lengths(mut self, encoding: Encoding<Integer>) -> Self {
        self.template_lengths = Some(encoding);
        self
    }

    pub fn set_mate_distances(mut self, encoding: Encoding<Integer>) -> Self {
        self.mate_distances = Some(encoding);
        self
    }

    pub fn set_tag_set_ids(mut self, encoding: Encoding<Integer>) -> Self {
        self.tag_set_ids = Some(encoding);
        self
    }

    pub fn set_feature_counts(mut self, encoding: Encoding<Integer>) -> Self {
        self.feature_counts = Some(encoding);
        self
    }

    pub fn set_feature_codes(mut self, encoding: Encoding<Byte>) -> Self {
        self.feature_codes = Some(encoding);
        self
    }

    pub fn set_feature_position_deltas(mut self, encoding: Encoding<Integer>) -> Self {
        self.feature_position_deltas = Some(encoding);
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

    pub fn set_insertion_bases(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.insertion_bases = Some(encoding);
        self
    }

    pub fn set_reference_skip_lengths(mut self, encoding: Encoding<Integer>) -> Self {
        self.reference_skip_lengths = Some(encoding);
        self
    }

    pub fn set_padding_lengths(mut self, encoding: Encoding<Integer>) -> Self {
        self.padding_lengths = Some(encoding);
        self
    }

    pub fn set_hard_clip_lengths(mut self, encoding: Encoding<Integer>) -> Self {
        self.hard_clip_lengths = Some(encoding);
        self
    }

    pub fn set_soft_clip_bases(mut self, encoding: Encoding<ByteArray>) -> Self {
        self.soft_clip_bases = Some(encoding);
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

    pub(crate) fn build(self) -> Result<DataSeriesEncodings, BuildError> {
        Ok(DataSeriesEncodings {
            bam_flags: self.bam_flags.ok_or(BuildError::MissingBamBitFlags)?,
            cram_flags: self.cram_flags.ok_or(BuildError::MissingCramBitFlags)?,
            reference_sequence_ids: self.reference_sequence_ids,
            read_lengths: self.read_lengths.ok_or(BuildError::MissingReadLengths)?,
            alignment_starts: self
                .alignment_starts
                .ok_or(BuildError::MissingInSeqPositions)?,
            read_group_ids: self.read_group_ids.ok_or(BuildError::MissingReadGroups)?,
            names: self.names,
            mate_flags: self.mate_flags,
            mate_reference_sequence_ids: self.mate_reference_sequence_ids,
            mate_alignment_starts: self.mate_alignment_starts,
            template_lengths: self.template_lengths,
            mate_distances: self.mate_distances,
            tag_set_ids: self.tag_set_ids.ok_or(BuildError::MissingTagIds)?,
            feature_counts: self.feature_counts,
            feature_codes: self.feature_codes,
            feature_position_deltas: self.feature_position_deltas,
            deletion_lengths: self.deletion_lengths,
            stretches_of_bases: self.stretches_of_bases,
            stretches_of_quality_scores: self.stretches_of_quality_scores,
            base_substitution_codes: self.base_substitution_codes,
            insertion_bases: self.insertion_bases,
            reference_skip_lengths: self.reference_skip_lengths,
            padding_lengths: self.padding_lengths,
            hard_clip_lengths: self.hard_clip_lengths,
            soft_clip_bases: self.soft_clip_bases,
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

        assert!(builder.bam_flags.is_none());
        assert!(builder.cram_flags.is_none());
        assert!(builder.reference_sequence_ids.is_none());
        assert!(builder.read_lengths.is_none());
        assert!(builder.alignment_starts.is_none());
        assert!(builder.read_group_ids.is_none());
        assert!(builder.names.is_none());
        assert!(builder.mate_flags.is_none());
        assert!(builder.mate_reference_sequence_ids.is_none());
        assert!(builder.mate_alignment_starts.is_none());
        assert!(builder.template_lengths.is_none());
        assert!(builder.mate_distances.is_none());
        assert!(builder.tag_set_ids.is_none());
        assert!(builder.feature_counts.is_none());
        assert!(builder.feature_codes.is_none());
        assert!(builder.feature_position_deltas.is_none());
        assert!(builder.deletion_lengths.is_none());
        assert!(builder.stretches_of_bases.is_none());
        assert!(builder.stretches_of_quality_scores.is_none());
        assert!(builder.base_substitution_codes.is_none());
        assert!(builder.insertion_bases.is_none());
        assert!(builder.reference_skip_lengths.is_none());
        assert!(builder.padding_lengths.is_none());
        assert!(builder.hard_clip_lengths.is_none());
        assert!(builder.soft_clip_bases.is_none());
        assert!(builder.mapping_qualities.is_none());
        assert!(builder.bases.is_none());
        assert!(builder.quality_scores.is_none());
    }
}
