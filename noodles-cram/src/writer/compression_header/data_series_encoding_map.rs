use std::io::{self, Write};

use crate::{
    container::compression_header::{data_series_encoding_map::DataSeries, DataSeriesEncodingMap},
    num::{write_itf8, Itf8},
    writer::encoding::write_encoding,
};

pub fn write_data_series_encoding_map<W>(
    writer: &mut W,
    data_series_encoding_map: &DataSeriesEncodingMap,
) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    let map_len = count_data_series_encodings(data_series_encoding_map);
    write_itf8(&mut buf, map_len)?;

    write_encodings(&mut buf, data_series_encoding_map)?;

    let data_len = buf.len() as Itf8;
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

fn count_data_series_encodings(data_series_encoding_map: &DataSeriesEncodingMap) -> Itf8 {
    // BAM bit flags, CRAM bit flags, read lengths, in-seq positions, read groups, tag IDs
    let mut n = 6;

    if data_series_encoding_map.reference_id_encoding().is_some() {
        n += 1;
    }

    if data_series_encoding_map.read_names_encoding().is_some() {
        n += 1;
    }

    if data_series_encoding_map
        .next_mate_bit_flags_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .next_fragment_reference_sequence_id_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .next_mate_alignment_start_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map.template_size_encoding().is_some() {
        n += 1;
    }

    if data_series_encoding_map
        .distance_to_next_fragment_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .number_of_read_features_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .read_features_codes_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .in_read_positions_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .deletion_lengths_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .stretches_of_bases_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .stretches_of_quality_scores_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map
        .base_substitution_codes_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map.insertion_encoding().is_some() {
        n += 1;
    }

    if data_series_encoding_map
        .reference_skip_length_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map.padding_encoding().is_some() {
        n += 1;
    }

    if data_series_encoding_map.hard_clip_encoding().is_some() {
        n += 1;
    }

    if data_series_encoding_map.soft_clip_encoding().is_some() {
        n += 1;
    }

    if data_series_encoding_map
        .mapping_qualities_encoding()
        .is_some()
    {
        n += 1;
    }

    if data_series_encoding_map.bases_encoding().is_some() {
        n += 1;
    }

    if data_series_encoding_map.quality_scores_encoding().is_some() {
        n += 1;
    }

    n
}

fn write_encodings<W>(
    writer: &mut W,
    data_series_encoding_map: &DataSeriesEncodingMap,
) -> io::Result<()>
where
    W: Write,
{
    write_key(writer, DataSeries::BamBitFlags)?;
    write_encoding(writer, data_series_encoding_map.bam_bit_flags_encoding())?;

    write_key(writer, DataSeries::CramBitFlags)?;
    write_encoding(writer, data_series_encoding_map.cram_bit_flags_encoding())?;

    if let Some(encoding) = data_series_encoding_map.reference_id_encoding() {
        write_key(writer, DataSeries::ReferenceId)?;
        write_encoding(writer, encoding)?;
    }

    write_key(writer, DataSeries::ReadLengths)?;
    write_encoding(writer, data_series_encoding_map.read_lengths_encoding())?;

    write_key(writer, DataSeries::InSeqPositions)?;
    write_encoding(writer, data_series_encoding_map.in_seq_positions_encoding())?;

    write_key(writer, DataSeries::ReadGroups)?;
    write_encoding(writer, data_series_encoding_map.read_groups_encoding())?;

    if let Some(encoding) = data_series_encoding_map.read_names_encoding() {
        write_key(writer, DataSeries::ReadNames)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.next_mate_bit_flags_encoding() {
        write_key(writer, DataSeries::NextMateBitFlags)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.next_fragment_reference_sequence_id_encoding()
    {
        write_key(writer, DataSeries::NextFragmentReferenceSequenceId)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.next_mate_alignment_start_encoding() {
        write_key(writer, DataSeries::NextMateAlignmentStart)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.template_size_encoding() {
        write_key(writer, DataSeries::TemplateSize)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.distance_to_next_fragment_encoding() {
        write_key(writer, DataSeries::DistanceToNextFragment)?;
        write_encoding(writer, encoding)?;
    }

    write_key(writer, DataSeries::TagIds)?;
    write_encoding(writer, data_series_encoding_map.tag_ids_encoding())?;

    if let Some(encoding) = data_series_encoding_map.number_of_read_features_encoding() {
        write_key(writer, DataSeries::NumberOfReadFeatures)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.read_features_codes_encoding() {
        write_key(writer, DataSeries::ReadFeaturesCodes)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.in_read_positions_encoding() {
        write_key(writer, DataSeries::InReadPositions)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.deletion_lengths_encoding() {
        write_key(writer, DataSeries::DeletionLengths)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.stretches_of_bases_encoding() {
        write_key(writer, DataSeries::StretchesOfBases)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.stretches_of_quality_scores_encoding() {
        write_key(writer, DataSeries::StretchesOfQualityScores)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.base_substitution_codes_encoding() {
        write_key(writer, DataSeries::BaseSubstitutionCodes)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.insertion_encoding() {
        write_key(writer, DataSeries::Insertion)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.reference_skip_length_encoding() {
        write_key(writer, DataSeries::ReferenceSkipLength)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.padding_encoding() {
        write_key(writer, DataSeries::Padding)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.hard_clip_encoding() {
        write_key(writer, DataSeries::HardClip)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.soft_clip_encoding() {
        write_key(writer, DataSeries::SoftClip)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.mapping_qualities_encoding() {
        write_key(writer, DataSeries::MappingQualities)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.bases_encoding() {
        write_key(writer, DataSeries::Bases)?;
        write_encoding(writer, encoding)?;
    }

    if let Some(encoding) = data_series_encoding_map.quality_scores_encoding() {
        write_key(writer, DataSeries::QualityScores)?;
        write_encoding(writer, encoding)?;
    }

    Ok(())
}
