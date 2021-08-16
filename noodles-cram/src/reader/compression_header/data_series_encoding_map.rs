use std::{
    convert::TryFrom,
    io::{self, Read},
};

use crate::{
    container::compression_header::{data_series_encoding_map::DataSeries, DataSeriesEncodingMap},
    num::read_itf8,
    reader::encoding::read_encoding,
};

pub fn read_data_series_encoding_map<R>(reader: &mut R) -> io::Result<DataSeriesEncodingMap>
where
    R: Read,
{
    let data_len = read_itf8(reader)?;
    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];
    let map_len = read_itf8(&mut buf_reader)?;

    let mut builder = DataSeriesEncodingMap::builder();
    let mut key_buf = [0; 2];

    for _ in 0..map_len {
        buf_reader.read_exact(&mut key_buf)?;

        let key = DataSeries::try_from(&key_buf[..])
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let encoding = read_encoding(&mut buf_reader)?;

        builder = match key {
            DataSeries::BamBitFlags => builder.set_bam_bit_flags_encoding(encoding),
            DataSeries::CramBitFlags => builder.set_cram_bit_flags_encoding(encoding),
            DataSeries::ReferenceId => builder.set_reference_id_encoding(encoding),
            DataSeries::ReadLengths => builder.set_read_lengths_encoding(encoding),
            DataSeries::InSeqPositions => builder.set_in_seq_positions_encoding(encoding),
            DataSeries::ReadGroups => builder.set_read_groups_encoding(encoding),
            DataSeries::ReadNames => builder.set_read_names_encoding(encoding),
            DataSeries::NextMateBitFlags => builder.set_next_mate_bit_flags_encoding(encoding),
            DataSeries::NextFragmentReferenceSequenceId => {
                builder.set_next_fragment_reference_sequence_id_encoding(encoding)
            }
            DataSeries::NextMateAlignmentStart => {
                builder.set_next_mate_alignment_start_encoding(encoding)
            }
            DataSeries::TemplateSize => builder.set_template_size_encoding(encoding),
            DataSeries::DistanceToNextFragment => {
                builder.set_distance_to_next_fragment_encoding(encoding)
            }
            DataSeries::TagIds => builder.set_tag_ids_encoding(encoding),
            DataSeries::NumberOfReadFeatures => {
                builder.set_number_of_read_features_encoding(encoding)
            }
            DataSeries::ReadFeaturesCodes => builder.set_read_features_codes_encoding(encoding),
            DataSeries::InReadPositions => builder.set_in_read_positions_encoding(encoding),
            DataSeries::DeletionLengths => builder.set_deletion_lengths_encoding(encoding),
            DataSeries::StretchesOfBases => builder.set_stretches_of_bases_encoding(encoding),
            DataSeries::StretchesOfQualityScores => {
                builder.set_stretches_of_quality_scores_encoding(encoding)
            }
            DataSeries::BaseSubstitutionCodes => {
                builder.set_base_substitution_codes_encoding(encoding)
            }
            DataSeries::Insertion => builder.set_insertion_encoding(encoding),
            DataSeries::ReferenceSkipLength => builder.set_reference_skip_length_encoding(encoding),
            DataSeries::Padding => builder.set_padding_encoding(encoding),
            DataSeries::HardClip => builder.set_hard_clip_encoding(encoding),
            DataSeries::SoftClip => builder.set_soft_clip_encoding(encoding),
            DataSeries::MappingQualities => builder.set_mapping_qualities_encoding(encoding),
            DataSeries::Bases => builder.set_bases_encoding(encoding),
            DataSeries::QualityScores => builder.set_quality_scores_encoding(encoding),
        }
    }

    Ok(builder.build())
}

#[cfg(test)]
mod tests {
    use crate::writer::compression_header::data_series_encoding_map::write_data_series_encoding_map;

    use super::*;

    fn build_data(data_series_encoding_map: &DataSeriesEncodingMap) -> io::Result<Vec<u8>> {
        let mut buf = Vec::new();
        write_data_series_encoding_map(&mut buf, data_series_encoding_map)?;
        Ok(buf)
    }

    #[test]
    fn test_read_data_series_encoding_map() -> io::Result<()> {
        let expected = DataSeriesEncodingMap::default();

        let data = build_data(&expected)?;
        let mut reader = &data[..];
        let actual = read_data_series_encoding_map(&mut reader)?;

        assert_eq!(
            actual.bam_bit_flags_encoding(),
            expected.bam_bit_flags_encoding()
        );

        assert_eq!(
            actual.cram_bit_flags_encoding(),
            expected.cram_bit_flags_encoding()
        );

        assert_eq!(
            actual.reference_id_encoding(),
            expected.reference_id_encoding()
        );

        assert_eq!(
            actual.read_lengths_encoding(),
            expected.read_lengths_encoding()
        );

        assert_eq!(
            actual.in_seq_positions_encoding(),
            expected.in_seq_positions_encoding()
        );

        assert_eq!(
            actual.read_groups_encoding(),
            expected.read_groups_encoding()
        );

        assert_eq!(actual.read_names_encoding(), expected.read_names_encoding());

        assert_eq!(
            actual.next_mate_bit_flags_encoding(),
            expected.next_mate_bit_flags_encoding()
        );

        assert_eq!(
            actual.next_fragment_reference_sequence_id_encoding(),
            expected.next_fragment_reference_sequence_id_encoding()
        );

        assert_eq!(
            actual.next_mate_alignment_start_encoding(),
            expected.next_mate_alignment_start_encoding()
        );

        assert_eq!(
            actual.template_size_encoding(),
            expected.template_size_encoding()
        );

        assert_eq!(
            actual.distance_to_next_fragment_encoding(),
            expected.distance_to_next_fragment_encoding()
        );

        assert_eq!(actual.tag_ids_encoding(), expected.tag_ids_encoding());

        assert_eq!(
            actual.number_of_read_features_encoding(),
            expected.number_of_read_features_encoding()
        );

        assert_eq!(
            actual.read_features_codes_encoding(),
            expected.read_features_codes_encoding()
        );

        assert_eq!(
            actual.in_read_positions_encoding(),
            expected.in_read_positions_encoding()
        );

        assert_eq!(
            actual.deletion_lengths_encoding(),
            expected.deletion_lengths_encoding()
        );

        assert_eq!(
            actual.stretches_of_bases_encoding(),
            expected.stretches_of_bases_encoding()
        );

        assert_eq!(
            actual.stretches_of_quality_scores_encoding(),
            expected.stretches_of_quality_scores_encoding()
        );

        assert_eq!(
            actual.base_substitution_codes_encoding(),
            expected.base_substitution_codes_encoding()
        );

        assert_eq!(actual.insertion_encoding(), expected.insertion_encoding());

        assert_eq!(
            actual.reference_skip_length_encoding(),
            expected.reference_skip_length_encoding()
        );

        assert_eq!(actual.padding_encoding(), expected.padding_encoding());
        assert_eq!(actual.hard_clip_encoding(), expected.hard_clip_encoding());
        assert_eq!(actual.soft_clip_encoding(), expected.soft_clip_encoding());

        assert_eq!(
            actual.mapping_qualities_encoding(),
            expected.mapping_qualities_encoding()
        );

        assert_eq!(actual.bases_encoding(), expected.bases_encoding());

        assert_eq!(
            actual.quality_scores_encoding(),
            expected.quality_scores_encoding()
        );

        Ok(())
    }
}
