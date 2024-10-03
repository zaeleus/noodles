use std::io::{self, Write};

use super::{
    write_encoding_for_byte_array_codec, write_encoding_for_byte_codec,
    write_encoding_for_integer_codec,
};
use crate::{
    data_container::compression_header::{
        data_series_encoding_map::DataSeries, DataSeriesEncodingMap,
    },
    io::writer::num::write_itf8,
};

pub(crate) fn write_data_series_encoding_map<W>(
    writer: &mut W,
    data_series_encoding_map: &DataSeriesEncodingMap,
) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    let map_len = i32::try_from(data_series_encoding_map.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut buf, map_len)?;

    write_encodings(&mut buf, data_series_encoding_map)?;

    let data_len =
        i32::try_from(buf.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, data_len)?;

    writer.write_all(&buf)
}

fn write_key<W>(writer: &mut W, key: DataSeries) -> io::Result<()>
where
    W: Write,
{
    let data = <[u8; 2]>::from(key);
    writer.write_all(&data)
}

fn write_encodings<W>(
    writer: &mut W,
    data_series_encoding_map: &DataSeriesEncodingMap,
) -> io::Result<()>
where
    W: Write,
{
    write_key(writer, DataSeries::BamFlags)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.bam_flags())?;

    write_key(writer, DataSeries::CramFlags)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.cram_flags())?;

    if let Some(encoding) = data_series_encoding_map.reference_sequence_ids() {
        write_key(writer, DataSeries::ReferenceSequenceIds)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    write_key(writer, DataSeries::ReadLengths)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.read_lengths())?;

    write_key(writer, DataSeries::AlignmentStarts)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.alignment_starts())?;

    write_key(writer, DataSeries::ReadGroupIds)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.read_group_ids())?;

    if let Some(encoding) = data_series_encoding_map.names() {
        write_key(writer, DataSeries::Names)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.mate_flags() {
        write_key(writer, DataSeries::MateFlags)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.mate_reference_sequence_ids() {
        write_key(writer, DataSeries::MateReferenceSequenceId)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.mate_alignment_starts() {
        write_key(writer, DataSeries::MateAlignmentStart)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.template_lengths() {
        write_key(writer, DataSeries::TemplateLengths)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.mate_distances() {
        write_key(writer, DataSeries::MateDistances)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    write_key(writer, DataSeries::TagSetIds)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.tag_set_ids())?;

    if let Some(encoding) = data_series_encoding_map.feature_counts() {
        write_key(writer, DataSeries::FeatureCounts)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.feature_codes() {
        write_key(writer, DataSeries::FeatureCodes)?;
        write_encoding_for_byte_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.feature_position_deltas() {
        write_key(writer, DataSeries::FeaturePositionDeltas)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.deletion_lengths() {
        write_key(writer, DataSeries::DeletionLengths)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.stretches_of_bases() {
        write_key(writer, DataSeries::StretchesOfBases)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.stretches_of_quality_scores() {
        write_key(writer, DataSeries::StretchesOfQualityScores)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.base_substitution_codes() {
        write_key(writer, DataSeries::BaseSubstitutionCodes)?;
        write_encoding_for_byte_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.insertion_bases() {
        write_key(writer, DataSeries::InsertionBases)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.reference_skip_lengths() {
        write_key(writer, DataSeries::ReferenceSkipLengths)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.padding_lengths() {
        write_key(writer, DataSeries::PaddingLengths)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.hard_clip_lengths() {
        write_key(writer, DataSeries::HardClipLengths)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.soft_clip_bases() {
        write_key(writer, DataSeries::SoftClipBases)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.mapping_qualities() {
        write_key(writer, DataSeries::MappingQualities)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.bases() {
        write_key(writer, DataSeries::Bases)?;
        write_encoding_for_byte_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.quality_scores() {
        write_key(writer, DataSeries::QualityScores)?;
        write_encoding_for_byte_codec(writer, encoding)?;
    }

    Ok(())
}
