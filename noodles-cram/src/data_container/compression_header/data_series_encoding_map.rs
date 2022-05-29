mod builder;
pub mod data_series;

pub use self::{builder::Builder, data_series::DataSeries};

use super::{
    encoding::codec::{Byte, ByteArray, Integer},
    Encoding,
};

/// A container compression header data series encoding map.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DataSeriesEncodingMap {
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
    pub fn builder() -> Builder {
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
            bam_bit_flags_encoding: Encoding::new(Integer::External(1)),
            cram_bit_flags_encoding: Encoding::new(Integer::External(2)),
            reference_id_encoding: Some(Encoding::new(Integer::External(3))),
            read_lengths_encoding: Encoding::new(Integer::External(4)),
            in_seq_positions_encoding: Encoding::new(Integer::External(5)),
            read_groups_encoding: Encoding::new(Integer::External(6)),
            read_names_encoding: Some(Encoding::new(ByteArray::ByteArrayStop(0x00, 7))),
            next_mate_bit_flags_encoding: Some(Encoding::new(Integer::External(8))),
            next_fragment_reference_sequence_id_encoding: Some(Encoding::new(Integer::External(9))),
            next_mate_alignment_start_encoding: Some(Encoding::new(Integer::External(10))),
            template_size_encoding: Some(Encoding::new(Integer::External(11))),
            distance_to_next_fragment_encoding: Some(Encoding::new(Integer::External(12))),
            tag_ids_encoding: Encoding::new(Integer::External(13)),
            number_of_read_features_encoding: Some(Encoding::new(Integer::External(14))),
            read_features_codes_encoding: Some(Encoding::new(Byte::External(15))),
            in_read_positions_encoding: Some(Encoding::new(Integer::External(16))),
            deletion_lengths_encoding: Some(Encoding::new(Integer::External(17))),
            stretches_of_bases_encoding: Some(Encoding::new(ByteArray::ByteArrayStop(0x00, 18))),
            stretches_of_quality_scores_encoding: Some(Encoding::new(ByteArray::ByteArrayLen(
                Encoding::new(Integer::External(19)),
                Encoding::new(Byte::External(19)),
            ))),
            base_substitution_codes_encoding: Some(Encoding::new(Byte::External(20))),
            insertion_encoding: Some(Encoding::new(ByteArray::ByteArrayStop(0x00, 21))),
            reference_skip_length_encoding: Some(Encoding::new(Integer::External(22))),
            padding_encoding: Some(Encoding::new(Integer::External(23))),
            hard_clip_encoding: Some(Encoding::new(Integer::External(24))),
            soft_clip_encoding: Some(Encoding::new(ByteArray::ByteArrayStop(0x00, 25))),
            mapping_qualities_encoding: Some(Encoding::new(Integer::External(26))),
            bases_encoding: Some(Encoding::new(Byte::External(27))),
            quality_scores_encoding: Some(Encoding::new(Byte::External(28))),
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
            .set_bam_bit_flags_encoding(Encoding::new(Integer::External(1)))
            .set_cram_bit_flags_encoding(Encoding::new(Integer::External(2)))
            .set_read_lengths_encoding(Encoding::new(Integer::External(4)))
            .set_in_seq_positions_encoding(Encoding::new(Integer::External(5)))
            .set_read_groups_encoding(Encoding::new(Integer::External(6)))
            .set_tag_ids_encoding(Encoding::new(Integer::External(13)))
            .build()?;

        assert_eq!(map.len(), 6);

        Ok(())
    }
}
