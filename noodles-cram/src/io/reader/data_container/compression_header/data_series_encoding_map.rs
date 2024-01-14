use std::io;

use bytes::{Buf, Bytes};

use super::{
    get_encoding_for_byte_array_codec, get_encoding_for_byte_codec, get_encoding_for_integer_codec,
};
use crate::{
    data_container::compression_header::{
        data_series_encoding_map::DataSeries, DataSeriesEncodingMap,
    },
    io::reader::num::get_itf8,
};

pub(super) fn get_data_series_encoding_map(src: &mut Bytes) -> io::Result<DataSeriesEncodingMap> {
    let data_len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < data_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut buf = src.split_to(data_len);

    let map_len = get_itf8(&mut buf)?;

    let mut builder = DataSeriesEncodingMap::builder();

    for _ in 0..map_len {
        let key = get_key(&mut buf)?;

        builder = match key {
            DataSeries::BamBitFlags => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_bam_bit_flags_encoding(encoding)
            }
            DataSeries::CramBitFlags => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_cram_bit_flags_encoding(encoding)
            }
            DataSeries::ReferenceId => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_reference_id_encoding(encoding)
            }
            DataSeries::ReadLengths => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_read_lengths_encoding(encoding)
            }
            DataSeries::InSeqPositions => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_in_seq_positions_encoding(encoding)
            }
            DataSeries::ReadGroups => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_read_groups_encoding(encoding)
            }
            DataSeries::ReadNames => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_read_names_encoding(encoding)
            }
            DataSeries::NextMateBitFlags => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_next_mate_bit_flags_encoding(encoding)
            }
            DataSeries::NextFragmentReferenceSequenceId => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_next_fragment_reference_sequence_id_encoding(encoding)
            }
            DataSeries::NextMateAlignmentStart => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_next_mate_alignment_start_encoding(encoding)
            }
            DataSeries::TemplateSize => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_template_size_encoding(encoding)
            }
            DataSeries::DistanceToNextFragment => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_distance_to_next_fragment_encoding(encoding)
            }
            DataSeries::TagIds => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_tag_ids_encoding(encoding)
            }
            DataSeries::NumberOfReadFeatures => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_number_of_read_features_encoding(encoding)
            }
            DataSeries::ReadFeaturesCodes => {
                let encoding = get_encoding_for_byte_codec(&mut buf)?;
                builder.set_read_features_codes_encoding(encoding)
            }
            DataSeries::InReadPositions => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_in_read_positions_encoding(encoding)
            }
            DataSeries::DeletionLengths => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_deletion_lengths_encoding(encoding)
            }
            DataSeries::StretchesOfBases => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_stretches_of_bases_encoding(encoding)
            }
            DataSeries::StretchesOfQualityScores => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_stretches_of_quality_scores_encoding(encoding)
            }
            DataSeries::BaseSubstitutionCodes => {
                let encoding = get_encoding_for_byte_codec(&mut buf)?;
                builder.set_base_substitution_codes_encoding(encoding)
            }
            DataSeries::Insertion => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_insertion_encoding(encoding)
            }
            DataSeries::ReferenceSkipLength => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_reference_skip_length_encoding(encoding)
            }
            DataSeries::Padding => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_padding_encoding(encoding)
            }
            DataSeries::HardClip => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_hard_clip_encoding(encoding)
            }
            DataSeries::SoftClip => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_soft_clip_encoding(encoding)
            }
            DataSeries::MappingQualities => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_mapping_qualities_encoding(encoding)
            }
            DataSeries::Bases => {
                let encoding = get_encoding_for_byte_codec(&mut buf)?;
                builder.set_bases_encoding(encoding)
            }
            DataSeries::QualityScores => {
                let encoding = get_encoding_for_byte_codec(&mut buf)?;
                builder.set_quality_scores_encoding(encoding)
            }
            DataSeries::ReservedTc | DataSeries::ReservedTn => {
                get_encoding_for_integer_codec(&mut buf)?;
                builder
            }
        }
    }

    builder
        .build()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn get_key<B>(src: &mut B) -> io::Result<DataSeries>
where
    B: Buf,
{
    let mut buf = [0; 2];

    if src.remaining() < buf.len() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    src.copy_to_slice(&mut buf);

    DataSeries::try_from(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_data(data_series_encoding_map: &DataSeriesEncodingMap) -> io::Result<Bytes> {
        use crate::io::writer::data_container::compression_header::data_series_encoding_map::write_data_series_encoding_map;

        let mut buf = Vec::new();
        write_data_series_encoding_map(&mut buf, data_series_encoding_map)?;
        Ok(Bytes::from(buf))
    }

    #[test]
    fn test_get_data_series_encoding_map() -> io::Result<()> {
        let expected = DataSeriesEncodingMap::default();

        let mut data = build_data(&expected)?;
        let actual = get_data_series_encoding_map(&mut data)?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
