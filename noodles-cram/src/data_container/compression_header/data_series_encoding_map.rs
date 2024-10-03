//! CRAM data container compress header data series-encoding map.

mod builder;
pub(crate) mod data_series;

pub(crate) use self::builder::Builder;
pub use self::data_series::DataSeries;

use super::{
    encoding::codec::{Byte, ByteArray, Integer},
    Encoding,
};
use crate::container::block;

/// A container compression header data series encoding map.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct DataSeriesEncodingMap {
    bam_flags: Encoding<Integer>,
    cram_flags: Encoding<Integer>,
    reference_sequence_ids: Option<Encoding<Integer>>,
    read_lengths: Encoding<Integer>,
    alignment_starts: Encoding<Integer>,
    read_group_ids: Encoding<Integer>,
    names: Option<Encoding<ByteArray>>,
    mate_flags: Option<Encoding<Integer>>,
    mate_reference_sequence_ids: Option<Encoding<Integer>>,
    mate_alignment_starts: Option<Encoding<Integer>>,
    template_lengths: Option<Encoding<Integer>>,
    mate_distances: Option<Encoding<Integer>>,
    tag_set_ids: Encoding<Integer>,
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

impl DataSeriesEncodingMap {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    pub fn len(&self) -> usize {
        // BAM bit flags, CRAM bit flags, read lengths, in-seq positions, read groups, tag IDs
        let mut n = 6;

        if self.reference_sequence_ids().is_some() {
            n += 1;
        }

        if self.names().is_some() {
            n += 1;
        }

        if self.mate_flags().is_some() {
            n += 1;
        }

        if self.mate_reference_sequence_ids().is_some() {
            n += 1;
        }

        if self.mate_alignment_starts().is_some() {
            n += 1;
        }

        if self.template_lengths().is_some() {
            n += 1;
        }

        if self.mate_distances().is_some() {
            n += 1;
        }

        if self.feature_counts().is_some() {
            n += 1;
        }

        if self.feature_codes().is_some() {
            n += 1;
        }

        if self.feature_position_deltas().is_some() {
            n += 1;
        }

        if self.deletion_lengths().is_some() {
            n += 1;
        }

        if self.stretches_of_bases().is_some() {
            n += 1;
        }

        if self.stretches_of_quality_scores().is_some() {
            n += 1;
        }

        if self.base_substitution_codes().is_some() {
            n += 1;
        }

        if self.insertion_bases().is_some() {
            n += 1;
        }

        if self.reference_skip_lengths().is_some() {
            n += 1;
        }

        if self.padding_lengths().is_some() {
            n += 1;
        }

        if self.hard_clip_lengths().is_some() {
            n += 1;
        }

        if self.soft_clip_bases().is_some() {
            n += 1;
        }

        if self.mapping_qualities().is_some() {
            n += 1;
        }

        if self.bases().is_some() {
            n += 1;
        }

        if self.quality_scores().is_some() {
            n += 1;
        }

        n
    }

    pub fn bam_flags(&self) -> &Encoding<Integer> {
        &self.bam_flags
    }

    pub fn cram_flags(&self) -> &Encoding<Integer> {
        &self.cram_flags
    }

    pub fn reference_sequence_ids(&self) -> Option<&Encoding<Integer>> {
        self.reference_sequence_ids.as_ref()
    }

    pub fn read_lengths(&self) -> &Encoding<Integer> {
        &self.read_lengths
    }

    pub fn alignment_starts(&self) -> &Encoding<Integer> {
        &self.alignment_starts
    }

    pub fn read_group_ids(&self) -> &Encoding<Integer> {
        &self.read_group_ids
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

    pub fn tag_set_ids(&self) -> &Encoding<Integer> {
        &self.tag_set_ids
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
}

impl Default for DataSeriesEncodingMap {
    fn default() -> Self {
        Self {
            bam_flags: Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::BamFlags),
            }),
            cram_flags: Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::CramFlags),
            }),
            reference_sequence_ids: Some(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::ReferenceSequenceIds),
            })),
            read_lengths: Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::ReadLengths),
            }),
            alignment_starts: Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::AlignmentStarts),
            }),
            read_group_ids: Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(DataSeries::ReadGroupIds),
            }),
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
            tag_set_ids: Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(13),
            }),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_len() -> Result<(), builder::BuildError> {
        let map = DataSeriesEncodingMap::default();
        assert_eq!(map.len(), 28);

        let map = DataSeriesEncodingMap::builder()
            .set_bam_flags(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(1),
            }))
            .set_cram_flags(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(2),
            }))
            .set_read_lengths(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(4),
            }))
            .set_alignment_starts(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(5),
            }))
            .set_read_group_ids(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(6),
            }))
            .set_tag_set_ids(Encoding::new(Integer::External {
                block_content_id: block::ContentId::from(13),
            }))
            .build()?;

        assert_eq!(map.len(), 6);

        Ok(())
    }
}
