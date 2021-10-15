use std::convert::TryFrom;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use super::read_encoding;
use crate::{
    data_container::compression_header::{
        data_series_encoding_map::DataSeries, DataSeriesEncodingMap,
    },
    r#async::reader::num::read_itf8,
};

pub async fn read_data_series_encoding_map<R>(reader: &mut R) -> io::Result<DataSeriesEncodingMap>
where
    R: AsyncRead + Unpin,
{
    let data_len = read_itf8(reader).await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = vec![0; data_len];
    reader.read_exact(&mut buf).await?;

    let mut buf_reader = &buf[..];
    read_data_series_encoding_map_map(&mut buf_reader).await
}

pub async fn read_data_series_encoding_map_map<R>(
    reader: &mut R,
) -> io::Result<DataSeriesEncodingMap>
where
    R: AsyncRead + Unpin,
{
    let map_len = read_itf8(reader).await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut builder = DataSeriesEncodingMap::builder();
    let mut key_buf = [0; 2];

    for _ in 0..map_len {
        reader.read_exact(&mut key_buf).await?;

        let key = DataSeries::try_from(key_buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let encoding = read_encoding(reader).await?;

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
            DataSeries::ReservedTc | DataSeries::ReservedTn => builder,
        }
    }

    builder
        .build()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_data_series_encoding_map() -> io::Result<()> {
        use crate::data_container::compression_header::Encoding;

        let data = [
            0x1f, // data.len = 31
            0x06, // map.len = 6
            b'B', b'F', 0x01, 0x01, 0x01, // (BamBitFlags, External(1))
            b'C', b'F', 0x01, 0x01, 0x02, // (CramBitFlags, External(2))
            b'R', b'L', 0x01, 0x01, 0x03, // (ReadLengths, External(3))
            b'A', b'P', 0x01, 0x01, 0x04, // (InSeqPositions, External(4))
            b'R', b'G', 0x01, 0x01, 0x05, // (ReadGroups, External(5))
            b'T', b'L', 0x01, 0x01, 0x06, // (TagIds, External(6))
        ];

        let mut reader = &data[..];
        let actual = read_data_series_encoding_map(&mut reader).await?;

        let expected = DataSeriesEncodingMap::builder()
            .set_bam_bit_flags_encoding(Encoding::External(1))
            .set_cram_bit_flags_encoding(Encoding::External(2))
            .set_read_lengths_encoding(Encoding::External(3))
            .set_in_seq_positions_encoding(Encoding::External(4))
            .set_read_groups_encoding(Encoding::External(5))
            .set_tag_ids_encoding(Encoding::External(6))
            .build()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
