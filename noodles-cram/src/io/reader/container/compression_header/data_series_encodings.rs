use std::io;

use super::{
    encoding::consume_any_encoding, read_byte_array_encoding, read_byte_encoding,
    read_integer_encoding,
};
use crate::{
    container::compression_header::{DataSeriesEncodings, data_series_encodings::DataSeries},
    file_definition::Version,
    io::reader::collections::read_map,
};

pub(super) fn read_data_series_encodings(
    src: &mut &[u8],
    version: Version,
) -> io::Result<DataSeriesEncodings> {
    let (mut buf, len) = read_map(src, version)?;
    read_data_series_encodings_inner(&mut buf, len, version)
}

fn read_data_series_encodings_inner(
    src: &mut &[u8],
    len: usize,
    version: Version,
) -> io::Result<DataSeriesEncodings> {
    let mut map = DataSeriesEncodings::default();

    for _ in 0..len {
        match read_key(src)? {
            DataSeries::BamFlags => {
                map.bam_flags = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::CramFlags => {
                map.cram_flags = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::ReferenceSequenceIds => {
                map.reference_sequence_ids = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::ReadLengths => {
                map.read_lengths = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::AlignmentStarts => {
                map.alignment_starts = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::ReadGroupIds => {
                map.read_group_ids = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::Names => map.names = read_byte_array_encoding(src, version).map(Some)?,
            DataSeries::MateFlags => {
                map.mate_flags = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::MateReferenceSequenceIds => {
                map.mate_reference_sequence_ids = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::MateAlignmentStarts => {
                map.mate_alignment_starts = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::TemplateLengths => {
                map.template_lengths = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::MateDistances => {
                map.mate_distances = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::TagSetIds => {
                map.tag_set_ids = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::FeatureCounts => {
                map.feature_counts = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::FeatureCodes => {
                map.feature_codes = read_byte_encoding(src, version).map(Some)?
            }
            DataSeries::FeaturePositionDeltas => {
                map.feature_position_deltas = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::DeletionLengths => {
                map.deletion_lengths = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::StretchesOfBases => {
                map.stretches_of_bases = read_byte_array_encoding(src, version).map(Some)?
            }
            DataSeries::StretchesOfQualityScores => {
                map.stretches_of_quality_scores =
                    read_byte_array_encoding(src, version).map(Some)?
            }
            DataSeries::BaseSubstitutionCodes => {
                map.base_substitution_codes = read_byte_encoding(src, version).map(Some)?
            }
            DataSeries::InsertionBases => {
                map.insertion_bases = read_byte_array_encoding(src, version).map(Some)?
            }
            DataSeries::ReferenceSkipLengths => {
                map.reference_skip_lengths = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::PaddingLengths => {
                map.padding_lengths = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::HardClipLengths => {
                map.hard_clip_lengths = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::SoftClipBases => {
                map.soft_clip_bases = read_byte_array_encoding(src, version).map(Some)?
            }
            DataSeries::MappingQualities => {
                map.mapping_qualities = read_integer_encoding(src, version).map(Some)?
            }
            DataSeries::Bases => map.bases = read_byte_encoding(src, version).map(Some)?,
            DataSeries::QualityScores => {
                map.quality_scores = read_byte_encoding(src, version).map(Some)?
            }
            DataSeries::ReservedTc | DataSeries::ReservedTn => {
                // § 8.4.2 "Compression header block: Data series encodings" (2024-09-04): "TC and
                // TN are legacy data series from CRAM 1.0... [and] decoders must silently skip
                // these fields."
                consume_any_encoding(src, version)?;
            }
        }
    }

    Ok(map)
}

fn read_key(src: &mut &[u8]) -> io::Result<DataSeries> {
    let (buf, rest) = src
        .split_first_chunk()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    DataSeries::try_from(*buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_data(data_series_encodings: &DataSeriesEncodings) -> io::Result<Vec<u8>> {
        use crate::file_definition::Version;
        use crate::io::writer::container::compression_header::data_series_encodings::write_data_series_encodings;

        let mut buf = Vec::new();
        write_data_series_encodings(&mut buf, data_series_encodings, Version::default())?;
        Ok(buf)
    }

    #[test]
    fn test_read_data_series_encodings() -> io::Result<()> {
        let expected = DataSeriesEncodings::init(Version::default());

        let src = build_data(&expected)?;
        let actual = read_data_series_encodings(&mut &src[..], Version::V3_0)?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
