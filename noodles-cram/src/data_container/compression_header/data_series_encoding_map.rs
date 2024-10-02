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
    bam_bit_flags_encoding: Encoding<Integer>,
    cram_bit_flags_encoding: Encoding<Integer>,
    reference_id_encoding: Option<Encoding<Integer>>,
    read_lengths_encoding: Encoding<Integer>,
    in_seq_positions_encoding: Encoding<Integer>,
    read_groups_encoding: Encoding<Integer>,
    read_names_encoding: Option<Encoding<ByteArray>>,
    next_mate_bit_flags_encoding: Option<Encoding<Integer>>,
    next_fragment_reference_sequence_id_encoding: Option<Encoding<Integer>>,
    next_mate_alignment_start_encoding: Option<Encoding<Integer>>,
    template_size_encoding: Option<Encoding<Integer>>,
    distance_to_next_fragment_encoding: Option<Encoding<Integer>>,
    tag_ids_encoding: Encoding<Integer>,
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

impl DataSeriesEncodingMap {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    pub fn len(&self) -> usize {
        // BAM bit flags, CRAM bit flags, read lengths, in-seq positions, read groups, tag IDs
        let mut n = 6;

        if self.reference_id_encoding().is_some() {
            n += 1;
        }

        if self.read_names_encoding().is_some() {
            n += 1;
        }

        if self.next_mate_bit_flags_encoding().is_some() {
            n += 1;
        }

        if self
            .next_fragment_reference_sequence_id_encoding()
            .is_some()
        {
            n += 1;
        }

        if self.next_mate_alignment_start_encoding().is_some() {
            n += 1;
        }

        if self.template_size_encoding().is_some() {
            n += 1;
        }

        if self.distance_to_next_fragment_encoding().is_some() {
            n += 1;
        }

        if self.number_of_read_features_encoding().is_some() {
            n += 1;
        }

        if self.read_features_codes_encoding().is_some() {
            n += 1;
        }

        if self.in_read_positions_encoding().is_some() {
            n += 1;
        }

        if self.deletion_lengths_encoding().is_some() {
            n += 1;
        }

        if self.stretches_of_bases_encoding().is_some() {
            n += 1;
        }

        if self.stretches_of_quality_scores_encoding().is_some() {
            n += 1;
        }

        if self.base_substitution_codes_encoding().is_some() {
            n += 1;
        }

        if self.insertion_encoding().is_some() {
            n += 1;
        }

        if self.reference_skip_length_encoding().is_some() {
            n += 1;
        }

        if self.padding_encoding().is_some() {
            n += 1;
        }

        if self.hard_clip_encoding().is_some() {
            n += 1;
        }

        if self.soft_clip_encoding().is_some() {
            n += 1;
        }

        if self.mapping_qualities_encoding().is_some() {
            n += 1;
        }

        if self.bases_encoding().is_some() {
            n += 1;
        }

        if self.quality_scores_encoding().is_some() {
            n += 1;
        }

        n
    }

    pub fn bam_bit_flags_encoding(&self) -> &Encoding<Integer> {
        &self.bam_bit_flags_encoding
    }

    pub fn cram_bit_flags_encoding(&self) -> &Encoding<Integer> {
        &self.cram_bit_flags_encoding
    }

    pub fn reference_id_encoding(&self) -> Option<&Encoding<Integer>> {
        self.reference_id_encoding.as_ref()
    }

    pub fn read_lengths_encoding(&self) -> &Encoding<Integer> {
        &self.read_lengths_encoding
    }

    pub fn in_seq_positions_encoding(&self) -> &Encoding<Integer> {
        &self.in_seq_positions_encoding
    }

    pub fn read_groups_encoding(&self) -> &Encoding<Integer> {
        &self.read_groups_encoding
    }

    pub fn read_names_encoding(&self) -> Option<&Encoding<ByteArray>> {
        self.read_names_encoding.as_ref()
    }

    pub fn next_mate_bit_flags_encoding(&self) -> Option<&Encoding<Integer>> {
        self.next_mate_bit_flags_encoding.as_ref()
    }

    pub fn next_fragment_reference_sequence_id_encoding(&self) -> Option<&Encoding<Integer>> {
        self.next_fragment_reference_sequence_id_encoding.as_ref()
    }

    pub fn next_mate_alignment_start_encoding(&self) -> Option<&Encoding<Integer>> {
        self.next_mate_alignment_start_encoding.as_ref()
    }

    pub fn template_size_encoding(&self) -> Option<&Encoding<Integer>> {
        self.template_size_encoding.as_ref()
    }

    pub fn distance_to_next_fragment_encoding(&self) -> Option<&Encoding<Integer>> {
        self.distance_to_next_fragment_encoding.as_ref()
    }

    pub fn tag_ids_encoding(&self) -> &Encoding<Integer> {
        &self.tag_ids_encoding
    }

    pub fn number_of_read_features_encoding(&self) -> Option<&Encoding<Integer>> {
        self.number_of_read_features_encoding.as_ref()
    }

    pub fn read_features_codes_encoding(&self) -> Option<&Encoding<Byte>> {
        self.read_features_codes_encoding.as_ref()
    }

    pub fn in_read_positions_encoding(&self) -> Option<&Encoding<Integer>> {
        self.in_read_positions_encoding.as_ref()
    }

    pub fn deletion_lengths_encoding(&self) -> Option<&Encoding<Integer>> {
        self.deletion_lengths_encoding.as_ref()
    }

    pub fn stretches_of_bases_encoding(&self) -> Option<&Encoding<ByteArray>> {
        self.stretches_of_bases_encoding.as_ref()
    }

    pub fn stretches_of_quality_scores_encoding(&self) -> Option<&Encoding<ByteArray>> {
        self.stretches_of_quality_scores_encoding.as_ref()
    }

    pub fn base_substitution_codes_encoding(&self) -> Option<&Encoding<Byte>> {
        self.base_substitution_codes_encoding.as_ref()
    }

    pub fn insertion_encoding(&self) -> Option<&Encoding<ByteArray>> {
        self.insertion_encoding.as_ref()
    }

    pub fn reference_skip_length_encoding(&self) -> Option<&Encoding<Integer>> {
        self.reference_skip_length_encoding.as_ref()
    }

    pub fn padding_encoding(&self) -> Option<&Encoding<Integer>> {
        self.padding_encoding.as_ref()
    }

    pub fn hard_clip_encoding(&self) -> Option<&Encoding<Integer>> {
        self.hard_clip_encoding.as_ref()
    }

    pub fn soft_clip_encoding(&self) -> Option<&Encoding<ByteArray>> {
        self.soft_clip_encoding.as_ref()
    }

    pub fn mapping_qualities_encoding(&self) -> Option<&Encoding<Integer>> {
        self.mapping_qualities_encoding.as_ref()
    }

    pub fn bases_encoding(&self) -> Option<&Encoding<Byte>> {
        self.bases_encoding.as_ref()
    }

    pub fn quality_scores_encoding(&self) -> Option<&Encoding<Byte>> {
        self.quality_scores_encoding.as_ref()
    }
}

impl Default for DataSeriesEncodingMap {
    fn default() -> Self {
        Self {
            bam_bit_flags_encoding: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::BamBitFlags,
            ))),
            cram_bit_flags_encoding: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::CramBitFlags,
            ))),
            reference_id_encoding: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::ReferenceId,
            )))),
            read_lengths_encoding: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::ReadLengths,
            ))),
            in_seq_positions_encoding: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::InSeqPositions,
            ))),
            read_groups_encoding: Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::ReadGroups,
            ))),
            read_names_encoding: Some(Encoding::new(ByteArray::ByteArrayStop(
                0x00,
                block::ContentId::from(DataSeries::ReadNames),
            ))),
            next_mate_bit_flags_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::NextMateBitFlags),
            ))),
            next_fragment_reference_sequence_id_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::NextFragmentReferenceSequenceId),
            ))),
            next_mate_alignment_start_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::NextMateAlignmentStart),
            ))),
            template_size_encoding: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::TemplateSize,
            )))),
            distance_to_next_fragment_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::DistanceToNextFragment),
            ))),
            tag_ids_encoding: Encoding::new(Integer::External(block::ContentId::from(13))),
            number_of_read_features_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::NumberOfReadFeatures),
            ))),
            read_features_codes_encoding: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::ReadFeaturesCodes),
            })),
            in_read_positions_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::InReadPositions),
            ))),
            deletion_lengths_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::DeletionLengths),
            ))),
            stretches_of_bases_encoding: Some(Encoding::new(ByteArray::ByteArrayStop(
                0x00,
                block::ContentId::from(DataSeries::StretchesOfBases),
            ))),
            stretches_of_quality_scores_encoding: Some(Encoding::new(ByteArray::ByteArrayLen(
                Encoding::new(Integer::External(block::ContentId::from(
                    DataSeries::StretchesOfQualityScores,
                ))),
                Encoding::new(Byte::External {
                    block_content_id: block::ContentId::from(DataSeries::StretchesOfQualityScores),
                }),
            ))),
            base_substitution_codes_encoding: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::BaseSubstitutionCodes),
            })),
            insertion_encoding: Some(Encoding::new(ByteArray::ByteArrayStop(
                0x00,
                block::ContentId::from(DataSeries::Insertion),
            ))),
            reference_skip_length_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::ReferenceSkipLength),
            ))),
            padding_encoding: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::Padding,
            )))),
            hard_clip_encoding: Some(Encoding::new(Integer::External(block::ContentId::from(
                DataSeries::HardClip,
            )))),
            soft_clip_encoding: Some(Encoding::new(ByteArray::ByteArrayStop(
                0x00,
                block::ContentId::from(DataSeries::SoftClip),
            ))),
            mapping_qualities_encoding: Some(Encoding::new(Integer::External(
                block::ContentId::from(DataSeries::MappingQualities),
            ))),
            bases_encoding: Some(Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::Bases),
            })),
            quality_scores_encoding: Some(Encoding::new(Byte::External {
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
            .set_bam_bit_flags_encoding(Encoding::new(Integer::External(block::ContentId::from(1))))
            .set_cram_bit_flags_encoding(Encoding::new(Integer::External(block::ContentId::from(
                2,
            ))))
            .set_read_lengths_encoding(Encoding::new(Integer::External(block::ContentId::from(4))))
            .set_in_seq_positions_encoding(Encoding::new(Integer::External(
                block::ContentId::from(5),
            )))
            .set_read_groups_encoding(Encoding::new(Integer::External(block::ContentId::from(6))))
            .set_tag_ids_encoding(Encoding::new(Integer::External(block::ContentId::from(13))))
            .build()?;

        assert_eq!(map.len(), 6);

        Ok(())
    }
}
