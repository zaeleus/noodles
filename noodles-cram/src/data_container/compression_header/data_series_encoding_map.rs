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
    bam_bit_flags: Encoding<Integer>,
    cram_bit_flags: Encoding<Integer>,
    reference_id: Option<Encoding<Integer>>,
    read_lengths: Encoding<Integer>,
    in_seq_positions: Encoding<Integer>,
    read_groups: Encoding<Integer>,
    read_names: Option<Encoding<ByteArray>>,
    next_mate_bit_flags: Option<Encoding<Integer>>,
    next_fragment_reference_sequence_id: Option<Encoding<Integer>>,
    next_mate_alignment_start: Option<Encoding<Integer>>,
    template_size: Option<Encoding<Integer>>,
    distance_to_next_fragment: Option<Encoding<Integer>>,
    tag_ids: Encoding<Integer>,
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

impl DataSeriesEncodingMap {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    pub fn len(&self) -> usize {
        // BAM bit flags, CRAM bit flags, read lengths, in-seq positions, read groups, tag IDs
        let mut n = 6;

        if self.reference_id().is_some() {
            n += 1;
        }

        if self.read_names().is_some() {
            n += 1;
        }

        if self.next_mate_bit_flags().is_some() {
            n += 1;
        }

        if self.next_fragment_reference_sequence_id().is_some() {
            n += 1;
        }

        if self.next_mate_alignment_start().is_some() {
            n += 1;
        }

        if self.template_size().is_some() {
            n += 1;
        }

        if self.distance_to_next_fragment().is_some() {
            n += 1;
        }

        if self.number_of_read_features().is_some() {
            n += 1;
        }

        if self.read_features_codes().is_some() {
            n += 1;
        }

        if self.in_read_positions().is_some() {
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

        if self.insertion().is_some() {
            n += 1;
        }

        if self.reference_skip_length().is_some() {
            n += 1;
        }

        if self.padding().is_some() {
            n += 1;
        }

        if self.hard_clip().is_some() {
            n += 1;
        }

        if self.soft_clip().is_some() {
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

    pub fn bam_bit_flags(&self) -> &Encoding<Integer> {
        &self.bam_bit_flags
    }

    pub fn cram_bit_flags(&self) -> &Encoding<Integer> {
        &self.cram_bit_flags
    }

    pub fn reference_id(&self) -> Option<&Encoding<Integer>> {
        self.reference_id.as_ref()
    }

    pub fn read_lengths(&self) -> &Encoding<Integer> {
        &self.read_lengths
    }

    pub fn in_seq_positions(&self) -> &Encoding<Integer> {
        &self.in_seq_positions
    }

    pub fn read_groups(&self) -> &Encoding<Integer> {
        &self.read_groups
    }

    pub fn read_names(&self) -> Option<&Encoding<ByteArray>> {
        self.read_names.as_ref()
    }

    pub fn next_mate_bit_flags(&self) -> Option<&Encoding<Integer>> {
        self.next_mate_bit_flags.as_ref()
    }

    pub fn next_fragment_reference_sequence_id(&self) -> Option<&Encoding<Integer>> {
        self.next_fragment_reference_sequence_id.as_ref()
    }

    pub fn next_mate_alignment_start(&self) -> Option<&Encoding<Integer>> {
        self.next_mate_alignment_start.as_ref()
    }

    pub fn template_size(&self) -> Option<&Encoding<Integer>> {
        self.template_size.as_ref()
    }

    pub fn distance_to_next_fragment(&self) -> Option<&Encoding<Integer>> {
        self.distance_to_next_fragment.as_ref()
    }

    pub fn tag_ids(&self) -> &Encoding<Integer> {
        &self.tag_ids
    }

    pub fn number_of_read_features(&self) -> Option<&Encoding<Integer>> {
        self.number_of_read_features.as_ref()
    }

    pub fn read_features_codes(&self) -> Option<&Encoding<Byte>> {
        self.read_features_codes.as_ref()
    }

    pub fn in_read_positions(&self) -> Option<&Encoding<Integer>> {
        self.in_read_positions.as_ref()
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

    pub fn insertion(&self) -> Option<&Encoding<ByteArray>> {
        self.insertion.as_ref()
    }

    pub fn reference_skip_length(&self) -> Option<&Encoding<Integer>> {
        self.reference_skip_length.as_ref()
    }

    pub fn padding(&self) -> Option<&Encoding<Integer>> {
        self.padding.as_ref()
    }

    pub fn hard_clip(&self) -> Option<&Encoding<Integer>> {
        self.hard_clip.as_ref()
    }

    pub fn soft_clip(&self) -> Option<&Encoding<ByteArray>> {
        self.soft_clip.as_ref()
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
            bam_bit_flags: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::BamBitFlags,
            ))),
            cram_bit_flags: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::CramBitFlags,
            ))),
            reference_id: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::ReferenceId,
            )))),
            read_lengths: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::ReadLengths,
            ))),
            in_seq_positions: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::InSeqPositions,
            ))),
            read_groups: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::ReadGroups,
            ))),
            read_names: Some(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: block::ContentId::from(DataSeries::ReadNames),
            })),
            next_mate_bit_flags: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::NextMateBitFlags,
            )))),
            next_fragment_reference_sequence_id: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::NextFragmentReferenceSequenceId),
            ))),
            next_mate_alignment_start: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::NextMateAlignmentStart),
            ))),
            template_size: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::TemplateSize,
            )))),
            distance_to_next_fragment: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::DistanceToNextFragment),
            ))),
            tag_ids: Encoding::new(Integer::External(block::ContentId::from(13))),
            number_of_read_features: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::NumberOfReadFeatures),
            ))),
            read_features_codes: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::ReadFeaturesCodes),
            })),
            in_read_positions: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::InReadPositions,
            )))),
            deletion_lengths: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::DeletionLengths,
            )))),
            stretches_of_bases: Some(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: block::ContentId::from(DataSeries::StretchesOfBases),
            })),
            stretches_of_quality_scores: Some(Encoding::new(ByteArray::ByteArrayLen {
                len_encoding: Encoding::new(Integer::External(block::ContentId::from(
                    DataSeries::StretchesOfQualityScores,
                ))),
                value_encoding: Encoding::new(Byte::External {
                    block_content_id: block::ContentId::from(DataSeries::StretchesOfQualityScores),
                }),
            })),
            base_substitution_codes: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::BaseSubstitutionCodes),
            })),
            insertion: Some(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: block::ContentId::from(DataSeries::Insertion),
            })),
            reference_skip_length: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::ReferenceSkipLength,
            )))),
            padding: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::Padding,
            )))),
            hard_clip: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::HardClip,
            )))),
            soft_clip: Some(Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: block::ContentId::from(DataSeries::SoftClip),
            })),
            mapping_qualities: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::MappingQualities,
            )))),
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
            .set_bam_bit_flags(Encoding::new(Integer::External(block::ContentId::from(1))))
            .set_cram_bit_flags(Encoding::new(Integer::External(block::ContentId::from(2))))
            .set_read_lengths(Encoding::new(Integer::External(block::ContentId::from(4))))
            .set_in_seq_positions(Encoding::new(Integer::External(block::ContentId::from(5))))
            .set_read_groups(Encoding::new(Integer::External(block::ContentId::from(6))))
            .set_tag_ids(Encoding::new(Integer::External(block::ContentId::from(13))))
            .build()?;

        assert_eq!(map.len(), 6);

        Ok(())
    }
}
