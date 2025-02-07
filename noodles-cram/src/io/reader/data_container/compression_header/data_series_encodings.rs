use std::io;

use bytes::{Buf, Bytes};

use super::{
    get_encoding_for_byte_array_codec, get_encoding_for_byte_codec, get_encoding_for_integer_codec,
};
use crate::{
    data_container::compression_header::{data_series_encodings::DataSeries, DataSeriesEncodings},
    io::reader::num::get_itf8,
};

pub(super) fn get_data_series_encodings(src: &mut Bytes) -> io::Result<DataSeriesEncodings> {
    let data_len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < data_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut buf = src.split_to(data_len);

    let map_len = get_itf8(&mut buf)?;

    let mut builder = DataSeriesEncodings::builder();

    for _ in 0..map_len {
        let key = get_key(&mut buf)?;

        builder = match key {
            DataSeries::BamFlags => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_bam_flags(encoding)
            }
            DataSeries::CramFlags => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_cram_flags(encoding)
            }
            DataSeries::ReferenceSequenceIds => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_reference_sequence_ids(encoding)
            }
            DataSeries::ReadLengths => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_read_lengths(encoding)
            }
            DataSeries::AlignmentStarts => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_alignment_starts(encoding)
            }
            DataSeries::ReadGroupIds => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_read_group_ids(encoding)
            }
            DataSeries::Names => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_names(encoding)
            }
            DataSeries::MateFlags => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_mate_flags(encoding)
            }
            DataSeries::MateReferenceSequenceId => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_mate_reference_sequence_ids(encoding)
            }
            DataSeries::MateAlignmentStart => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_mate_alignment_starts(encoding)
            }
            DataSeries::TemplateLengths => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_template_lengths(encoding)
            }
            DataSeries::MateDistances => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_mate_distances(encoding)
            }
            DataSeries::TagSetIds => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_tag_set_ids(encoding)
            }
            DataSeries::FeatureCounts => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_feature_counts(encoding)
            }
            DataSeries::FeatureCodes => {
                let encoding = get_encoding_for_byte_codec(&mut buf)?;
                builder.set_feature_codes(encoding)
            }
            DataSeries::FeaturePositionDeltas => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_feature_position_deltas(encoding)
            }
            DataSeries::DeletionLengths => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_deletion_lengths(encoding)
            }
            DataSeries::StretchesOfBases => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_stretches_of_bases(encoding)
            }
            DataSeries::StretchesOfQualityScores => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_stretches_of_quality_scores(encoding)
            }
            DataSeries::BaseSubstitutionCodes => {
                let encoding = get_encoding_for_byte_codec(&mut buf)?;
                builder.set_base_substitution_codes(encoding)
            }
            DataSeries::InsertionBases => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_insertion_bases(encoding)
            }
            DataSeries::ReferenceSkipLengths => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_reference_skip_lengths(encoding)
            }
            DataSeries::PaddingLengths => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_padding_lengths(encoding)
            }
            DataSeries::HardClipLengths => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_hard_clip_lengths(encoding)
            }
            DataSeries::SoftClipBases => {
                let encoding = get_encoding_for_byte_array_codec(&mut buf)?;
                builder.set_soft_clip_bases(encoding)
            }
            DataSeries::MappingQualities => {
                let encoding = get_encoding_for_integer_codec(&mut buf)?;
                builder.set_mapping_qualities(encoding)
            }
            DataSeries::Bases => {
                let encoding = get_encoding_for_byte_codec(&mut buf)?;
                builder.set_bases(encoding)
            }
            DataSeries::QualityScores => {
                let encoding = get_encoding_for_byte_codec(&mut buf)?;
                builder.set_quality_scores(encoding)
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

    fn build_data(data_series_encodings: &DataSeriesEncodings) -> io::Result<Bytes> {
        use crate::io::writer::data_container::compression_header::data_series_encodings::write_data_series_encodings;

        let mut buf = Vec::new();
        write_data_series_encodings(&mut buf, data_series_encodings)?;
        Ok(Bytes::from(buf))
    }

    #[test]
    fn test_get_data_series_encodings() -> io::Result<()> {
        let expected = DataSeriesEncodings::default();

        let mut data = build_data(&expected)?;
        let actual = get_data_series_encodings(&mut data)?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
