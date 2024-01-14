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
    write_key(writer, DataSeries::BamBitFlags)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.bam_bit_flags_encoding())?;

    write_key(writer, DataSeries::CramBitFlags)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.cram_bit_flags_encoding())?;

    if let Some(encoding) = data_series_encoding_map.reference_id_encoding() {
        write_key(writer, DataSeries::ReferenceId)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    write_key(writer, DataSeries::ReadLengths)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.read_lengths_encoding())?;

    write_key(writer, DataSeries::InSeqPositions)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.in_seq_positions_encoding())?;

    write_key(writer, DataSeries::ReadGroups)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.read_groups_encoding())?;

    if let Some(encoding) = data_series_encoding_map.read_names_encoding() {
        write_key(writer, DataSeries::ReadNames)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.next_mate_bit_flags_encoding() {
        write_key(writer, DataSeries::NextMateBitFlags)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.next_fragment_reference_sequence_id_encoding()
    {
        write_key(writer, DataSeries::NextFragmentReferenceSequenceId)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.next_mate_alignment_start_encoding() {
        write_key(writer, DataSeries::NextMateAlignmentStart)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.template_size_encoding() {
        write_key(writer, DataSeries::TemplateSize)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.distance_to_next_fragment_encoding() {
        write_key(writer, DataSeries::DistanceToNextFragment)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    write_key(writer, DataSeries::TagIds)?;
    write_encoding_for_integer_codec(writer, data_series_encoding_map.tag_ids_encoding())?;

    if let Some(encoding) = data_series_encoding_map.number_of_read_features_encoding() {
        write_key(writer, DataSeries::NumberOfReadFeatures)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.read_features_codes_encoding() {
        write_key(writer, DataSeries::ReadFeaturesCodes)?;
        write_encoding_for_byte_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.in_read_positions_encoding() {
        write_key(writer, DataSeries::InReadPositions)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.deletion_lengths_encoding() {
        write_key(writer, DataSeries::DeletionLengths)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.stretches_of_bases_encoding() {
        write_key(writer, DataSeries::StretchesOfBases)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.stretches_of_quality_scores_encoding() {
        write_key(writer, DataSeries::StretchesOfQualityScores)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.base_substitution_codes_encoding() {
        write_key(writer, DataSeries::BaseSubstitutionCodes)?;
        write_encoding_for_byte_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.insertion_encoding() {
        write_key(writer, DataSeries::Insertion)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.reference_skip_length_encoding() {
        write_key(writer, DataSeries::ReferenceSkipLength)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.padding_encoding() {
        write_key(writer, DataSeries::Padding)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.hard_clip_encoding() {
        write_key(writer, DataSeries::HardClip)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.soft_clip_encoding() {
        write_key(writer, DataSeries::SoftClip)?;
        write_encoding_for_byte_array_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.mapping_qualities_encoding() {
        write_key(writer, DataSeries::MappingQualities)?;
        write_encoding_for_integer_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.bases_encoding() {
        write_key(writer, DataSeries::Bases)?;
        write_encoding_for_byte_codec(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.quality_scores_encoding() {
        write_key(writer, DataSeries::QualityScores)?;
        write_encoding_for_byte_codec(writer, encoding)?;
    }

    Ok(())
}
