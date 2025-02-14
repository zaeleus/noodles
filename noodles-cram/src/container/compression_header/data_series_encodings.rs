//! CRAM container compression header data series encodings.

pub(crate) mod data_series;

pub use self::data_series::DataSeries;

use super::{
    encoding::codec::{Byte, ByteArray, Integer},
    Encoding,
};
use crate::container::block;

/// CRAM container compression header data series encodings.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct DataSeriesEncodings {
    pub(crate) bam_flags: Option<Encoding<Integer>>,
    pub(crate) cram_flags: Option<Encoding<Integer>>,
    pub(crate) reference_sequence_ids: Option<Encoding<Integer>>,
    pub(crate) read_lengths: Option<Encoding<Integer>>,
    pub(crate) alignment_starts: Option<Encoding<Integer>>,
    pub(crate) read_group_ids: Option<Encoding<Integer>>,
    pub(crate) names: Option<Encoding<ByteArray>>,
    pub(crate) mate_flags: Option<Encoding<Integer>>,
    pub(crate) mate_reference_sequence_ids: Option<Encoding<Integer>>,
    pub(crate) mate_alignment_starts: Option<Encoding<Integer>>,
    pub(crate) template_lengths: Option<Encoding<Integer>>,
    pub(crate) mate_distances: Option<Encoding<Integer>>,
    pub(crate) tag_set_ids: Option<Encoding<Integer>>,
    pub(crate) feature_counts: Option<Encoding<Integer>>,
    pub(crate) feature_codes: Option<Encoding<Byte>>,
    pub(crate) feature_position_deltas: Option<Encoding<Integer>>,
    pub(crate) deletion_lengths: Option<Encoding<Integer>>,
    pub(crate) stretches_of_bases: Option<Encoding<ByteArray>>,
    pub(crate) stretches_of_quality_scores: Option<Encoding<ByteArray>>,
    pub(crate) base_substitution_codes: Option<Encoding<Byte>>,
    pub(crate) insertion_bases: Option<Encoding<ByteArray>>,
    pub(crate) reference_skip_lengths: Option<Encoding<Integer>>,
    pub(crate) padding_lengths: Option<Encoding<Integer>>,
    pub(crate) hard_clip_lengths: Option<Encoding<Integer>>,
    pub(crate) soft_clip_bases: Option<Encoding<ByteArray>>,
    pub(crate) mapping_qualities: Option<Encoding<Integer>>,
    pub(crate) bases: Option<Encoding<Byte>>,
    pub(crate) quality_scores: Option<Encoding<Byte>>,
}

impl DataSeriesEncodings {
    pub fn bam_flags(&self) -> Option<&Encoding<Integer>> {
        self.bam_flags.as_ref()
    }

    pub fn cram_flags(&self) -> Option<&Encoding<Integer>> {
        self.cram_flags.as_ref()
    }

    pub fn reference_sequence_ids(&self) -> Option<&Encoding<Integer>> {
        self.reference_sequence_ids.as_ref()
    }

    pub fn read_lengths(&self) -> Option<&Encoding<Integer>> {
        self.read_lengths.as_ref()
    }

    pub fn alignment_starts(&self) -> Option<&Encoding<Integer>> {
        self.alignment_starts.as_ref()
    }

    pub fn read_group_ids(&self) -> Option<&Encoding<Integer>> {
        self.read_group_ids.as_ref()
    }

    pub fn names(&self) -> Option<&Encoding<ByteArray>> {
        self.names.as_ref()
    }

    pub fn mate_flags(&self) -> Option<&Encoding<Integer>> {
        self.mate_flags.as_ref()
    }

    pub fn mate_reference_sequence_ids(&self) -> Option<&Encoding<Integer>> {
        self.mate_reference_sequence_ids.as_ref()
    }

    pub fn mate_alignment_starts(&self) -> Option<&Encoding<Integer>> {
        self.mate_alignment_starts.as_ref()
    }

    pub fn template_lengths(&self) -> Option<&Encoding<Integer>> {
        self.template_lengths.as_ref()
    }

    pub fn mate_distances(&self) -> Option<&Encoding<Integer>> {
        self.mate_distances.as_ref()
    }

    pub fn tag_set_ids(&self) -> Option<&Encoding<Integer>> {
        self.tag_set_ids.as_ref()
    }

    pub fn feature_counts(&self) -> Option<&Encoding<Integer>> {
        self.feature_counts.as_ref()
    }

    pub fn feature_codes(&self) -> Option<&Encoding<Byte>> {
        self.feature_codes.as_ref()
    }

    pub fn feature_position_deltas(&self) -> Option<&Encoding<Integer>> {
        self.feature_position_deltas.as_ref()
    }

    pub fn deletion_lengths(&self) -> Option<&Encoding<Integer>> {
        self.deletion_lengths.as_ref()
    }

    pub fn stretches_of_bases(&self) -> Option<&Encoding<ByteArray>> {
        self.stretches_of_bases.as_ref()
    }

    pub fn stretches_of_quality_scores(&self) -> Option<&Encoding<ByteArray>> {
        self.stretches_of_quality_scores.as_ref()
    }

    pub fn base_substitution_codes(&self) -> Option<&Encoding<Byte>> {
        self.base_substitution_codes.as_ref()
    }

    pub fn insertion_bases(&self) -> Option<&Encoding<ByteArray>> {
        self.insertion_bases.as_ref()
    }

    pub fn reference_skip_lengths(&self) -> Option<&Encoding<Integer>> {
        self.reference_skip_lengths.as_ref()
    }

    pub fn padding_lengths(&self) -> Option<&Encoding<Integer>> {
        self.padding_lengths.as_ref()
    }

    pub fn hard_clip_lengths(&self) -> Option<&Encoding<Integer>> {
        self.hard_clip_lengths.as_ref()
    }

    pub fn soft_clip_bases(&self) -> Option<&Encoding<ByteArray>> {
        self.soft_clip_bases.as_ref()
    }

    pub fn mapping_qualities(&self) -> Option<&Encoding<Integer>> {
        self.mapping_qualities.as_ref()
    }

    pub fn bases(&self) -> Option<&Encoding<Byte>> {
        self.bases.as_ref()
    }

    pub fn quality_scores(&self) -> Option<&Encoding<Byte>> {
        self.quality_scores.as_ref()
    }

    pub fn init() -> Self {
        Self {
            bam_flags: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::BamFlags),
            })),
            cram_flags: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::CramFlags),
            })),
            reference_sequence_ids: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::ReferenceSequenceIds),
            })),
            read_lengths: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::ReadLengths),
            })),
            alignment_starts: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::AlignmentStarts),
            })),
            read_group_ids: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::ReadGroupIds),
            })),
            names: Some(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: block::ContentId::from(DataSeries::Names),
            })),
            mate_flags: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::MateFlags),
            })),
            mate_reference_sequence_ids: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::MateReferenceSequenceId),
            })),
            mate_alignment_starts: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::MateAlignmentStart),
            })),
            template_lengths: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::TemplateLengths),
            })),
            mate_distances: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::MateDistances),
            })),
            tag_set_ids: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::TagSetIds),
            })),
            feature_counts: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::FeatureCounts),
            })),
            feature_codes: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::FeatureCodes),
            })),
            feature_position_deltas: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::FeaturePositionDeltas),
            })),
            deletion_lengths: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::DeletionLengths),
            })),
            stretches_of_bases: Some(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: block::ContentId::from(DataSeries::StretchesOfBases),
            })),
            stretches_of_quality_scores: Some(Encoding::new(ByteArray::ByteArrayLen {
                len_encoding: Encoding::new(Integer::External {
                    block_content_id: block::ContentId::from(DataSeries::StretchesOfQualityScores),
                }),
                value_encoding: Encoding::new(Byte::External {
                    block_content_id: block::ContentId::from(DataSeries::StretchesOfQualityScores),
                }),
            })),
            base_substitution_codes: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::BaseSubstitutionCodes),
            })),
            insertion_bases: Some(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: block::ContentId::from(DataSeries::InsertionBases),
            })),
            reference_skip_lengths: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::ReferenceSkipLengths),
            })),
            padding_lengths: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::PaddingLengths),
            })),
            hard_clip_lengths: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::HardClipLengths),
            })),
            soft_clip_bases: Some(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: block::ContentId::from(DataSeries::SoftClipBases),
            })),
            mapping_qualities: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::MappingQualities),
            })),
            bases: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::Bases),
            })),
            quality_scores: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::QualityScores),
            })),
        }
    }
}
