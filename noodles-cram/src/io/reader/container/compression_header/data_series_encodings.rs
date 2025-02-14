use std::io;

use super::{read_byte_array_encoding, read_byte_encoding, read_integer_encoding};
use crate::{
    container::compression_header::{data_series_encodings::DataSeries, DataSeriesEncodings},
    io::reader::{collections::read_map, split_at_checked},
};

pub(super) fn read_data_series_encodings(src: &mut &[u8]) -> io::Result<DataSeriesEncodings> {
    let (mut buf, len) = read_map(src)?;
    read_data_series_encodings_inner(&mut buf, len)
}

fn read_data_series_encodings_inner(
    src: &mut &[u8],
    len: usize,
) -> io::Result<DataSeriesEncodings> {
    let mut builder = DataSeriesEncodings::builder();

    for _ in 0..len {
        let key = read_key(src)?;

        builder = match key {
            DataSeries::BamFlags => {
                let encoding = read_integer_encoding(src)?;
                builder.set_bam_flags(encoding)
            }
            DataSeries::CramFlags => {
                let encoding = read_integer_encoding(src)?;
                builder.set_cram_flags(encoding)
            }
            DataSeries::ReferenceSequenceIds => {
                let encoding = read_integer_encoding(src)?;
                builder.set_reference_sequence_ids(encoding)
            }
            DataSeries::ReadLengths => {
                let encoding = read_integer_encoding(src)?;
                builder.set_read_lengths(encoding)
            }
            DataSeries::AlignmentStarts => {
                let encoding = read_integer_encoding(src)?;
                builder.set_alignment_starts(encoding)
            }
            DataSeries::ReadGroupIds => {
                let encoding = read_integer_encoding(src)?;
                builder.set_read_group_ids(encoding)
            }
            DataSeries::Names => {
                let encoding = read_byte_array_encoding(src)?;
                builder.set_names(encoding)
            }
            DataSeries::MateFlags => {
                let encoding = read_integer_encoding(src)?;
                builder.set_mate_flags(encoding)
            }
            DataSeries::MateReferenceSequenceId => {
                let encoding = read_integer_encoding(src)?;
                builder.set_mate_reference_sequence_ids(encoding)
            }
            DataSeries::MateAlignmentStart => {
                let encoding = read_integer_encoding(src)?;
                builder.set_mate_alignment_starts(encoding)
            }
            DataSeries::TemplateLengths => {
                let encoding = read_integer_encoding(src)?;
                builder.set_template_lengths(encoding)
            }
            DataSeries::MateDistances => {
                let encoding = read_integer_encoding(src)?;
                builder.set_mate_distances(encoding)
            }
            DataSeries::TagSetIds => {
                let encoding = read_integer_encoding(src)?;
                builder.set_tag_set_ids(encoding)
            }
            DataSeries::FeatureCounts => {
                let encoding = read_integer_encoding(src)?;
                builder.set_feature_counts(encoding)
            }
            DataSeries::FeatureCodes => {
                let encoding = read_byte_encoding(src)?;
                builder.set_feature_codes(encoding)
            }
            DataSeries::FeaturePositionDeltas => {
                let encoding = read_integer_encoding(src)?;
                builder.set_feature_position_deltas(encoding)
            }
            DataSeries::DeletionLengths => {
                let encoding = read_integer_encoding(src)?;
                builder.set_deletion_lengths(encoding)
            }
            DataSeries::StretchesOfBases => {
                let encoding = read_byte_array_encoding(src)?;
                builder.set_stretches_of_bases(encoding)
            }
            DataSeries::StretchesOfQualityScores => {
                let encoding = read_byte_array_encoding(src)?;
                builder.set_stretches_of_quality_scores(encoding)
            }
            DataSeries::BaseSubstitutionCodes => {
                let encoding = read_byte_encoding(src)?;
                builder.set_base_substitution_codes(encoding)
            }
            DataSeries::InsertionBases => {
                let encoding = read_byte_array_encoding(src)?;
                builder.set_insertion_bases(encoding)
            }
            DataSeries::ReferenceSkipLengths => {
                let encoding = read_integer_encoding(src)?;
                builder.set_reference_skip_lengths(encoding)
            }
            DataSeries::PaddingLengths => {
                let encoding = read_integer_encoding(src)?;
                builder.set_padding_lengths(encoding)
            }
            DataSeries::HardClipLengths => {
                let encoding = read_integer_encoding(src)?;
                builder.set_hard_clip_lengths(encoding)
            }
            DataSeries::SoftClipBases => {
                let encoding = read_byte_array_encoding(src)?;
                builder.set_soft_clip_bases(encoding)
            }
            DataSeries::MappingQualities => {
                let encoding = read_integer_encoding(src)?;
                builder.set_mapping_qualities(encoding)
            }
            DataSeries::Bases => {
                let encoding = read_byte_encoding(src)?;
                builder.set_bases(encoding)
            }
            DataSeries::QualityScores => {
                let encoding = read_byte_encoding(src)?;
                builder.set_quality_scores(encoding)
            }
            DataSeries::ReservedTc | DataSeries::ReservedTn => {
                read_integer_encoding(src)?;
                builder
            }
        }
    }

    builder
        .build()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_key(src: &mut &[u8]) -> io::Result<DataSeries> {
    const SIZE: usize = 2;

    let (buf, rest) =
        split_at_checked(src, SIZE).ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    // SAFETY: `buf.len() == 2`.
    let key: [u8; SIZE] = buf.try_into().unwrap();

    DataSeries::try_from(key).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_data(data_series_encodings: &DataSeriesEncodings) -> io::Result<Vec<u8>> {
        use crate::io::writer::container::compression_header::data_series_encodings::write_data_series_encodings;

        let mut buf = Vec::new();
        write_data_series_encodings(&mut buf, data_series_encodings)?;
        Ok(buf)
    }

    #[test]
    fn test_read_data_series_encodings() -> io::Result<()> {
        let expected = DataSeriesEncodings::default();

        let src = build_data(&expected)?;
        let actual = read_data_series_encodings(&mut &src[..])?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
