use std::io::{self, Write};

use super::{write_byte_array_encoding, write_byte_encoding, write_integer_encoding};
use crate::{
    container::compression_header::{
        DataSeriesEncodings, Encoding,
        data_series_encodings::DataSeries,
        encoding::codec::{Byte, ByteArray, Integer},
    },
    file_definition::Version,
    io::writer::{collections::write_array, num::write_int},
};

pub(crate) fn write_data_series_encodings<W>(
    writer: &mut W,
    data_series_encodings: &DataSeriesEncodings,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let buf = encode(data_series_encodings, version)?;
    write_array(writer, version, &buf)
}

fn encode(data_series_encodings: &DataSeriesEncodings, version: Version) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    encode_inner(&mut buf, data_series_encodings, version)?;
    Ok(buf)
}

fn encode_inner<W>(
    writer: &mut W,
    data_series_encodings: &DataSeriesEncodings,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let len = data_series_encodings_len(data_series_encodings);
    let n = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, n)?;

    write_encodings(writer, data_series_encodings, version)?;

    Ok(())
}

fn data_series_encodings_len(encodings: &DataSeriesEncodings) -> usize {
    fn count(n: &mut usize, is_some: bool) {
        if is_some {
            *n += 1;
        }
    }

    let mut n = 0;

    count(&mut n, encodings.bam_flags().is_some());
    count(&mut n, encodings.cram_flags().is_some());
    count(&mut n, encodings.reference_sequence_ids().is_some());
    count(&mut n, encodings.read_lengths().is_some());
    count(&mut n, encodings.alignment_starts().is_some());
    count(&mut n, encodings.read_group_ids().is_some());
    count(&mut n, encodings.names().is_some());
    count(&mut n, encodings.mate_flags().is_some());
    count(&mut n, encodings.mate_reference_sequence_ids().is_some());
    count(&mut n, encodings.mate_alignment_starts().is_some());
    count(&mut n, encodings.template_lengths().is_some());
    count(&mut n, encodings.mate_distances().is_some());
    count(&mut n, encodings.tag_set_ids().is_some());
    count(&mut n, encodings.feature_counts().is_some());
    count(&mut n, encodings.feature_codes().is_some());
    count(&mut n, encodings.feature_position_deltas().is_some());
    count(&mut n, encodings.deletion_lengths().is_some());
    count(&mut n, encodings.stretches_of_bases().is_some());
    count(&mut n, encodings.stretches_of_quality_scores().is_some());
    count(&mut n, encodings.base_substitution_codes().is_some());
    count(&mut n, encodings.insertion_bases().is_some());
    count(&mut n, encodings.reference_skip_lengths().is_some());
    count(&mut n, encodings.padding_lengths().is_some());
    count(&mut n, encodings.hard_clip_lengths().is_some());
    count(&mut n, encodings.soft_clip_bases().is_some());
    count(&mut n, encodings.mapping_qualities().is_some());
    count(&mut n, encodings.bases().is_some());
    count(&mut n, encodings.quality_scores().is_some());

    n
}

fn write_key<W>(writer: &mut W, key: DataSeries) -> io::Result<()>
where
    W: Write,
{
    let data = <[u8; 2]>::from(key);
    writer.write_all(&data)
}

#[rustfmt::skip]
fn write_encodings<W>(writer: &mut W, encodings: &DataSeriesEncodings, version: Version) -> io::Result<()>
where
    W: Write,
{
    maybe_write_integer_encoding(writer, DataSeries::BamFlags, encodings.bam_flags(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::CramFlags, encodings.cram_flags(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::ReferenceSequenceIds, encodings.reference_sequence_ids(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::ReadLengths, encodings.read_lengths(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::AlignmentStarts, encodings.alignment_starts(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::ReadGroupIds, encodings.read_group_ids(), version)?;
    maybe_write_byte_array_encoding(writer, DataSeries::Names, encodings.names(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::MateFlags, encodings.mate_flags(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::MateReferenceSequenceIds, encodings.mate_reference_sequence_ids(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::MateAlignmentStarts, encodings.mate_alignment_starts(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::TemplateLengths, encodings.template_lengths(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::MateDistances, encodings.mate_distances(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::TagSetIds, encodings.tag_set_ids(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::FeatureCounts, encodings.feature_counts(), version)?;
    maybe_write_byte_encoding(writer, DataSeries::FeatureCodes, encodings.feature_codes(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::FeaturePositionDeltas, encodings.feature_position_deltas(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::DeletionLengths, encodings.deletion_lengths(), version)?;
    maybe_write_byte_array_encoding(writer, DataSeries::StretchesOfBases, encodings.stretches_of_bases(), version)?;
    maybe_write_byte_array_encoding(writer, DataSeries::StretchesOfQualityScores, encodings.stretches_of_quality_scores(), version)?;
    maybe_write_byte_encoding(writer, DataSeries::BaseSubstitutionCodes, encodings.base_substitution_codes(), version)?;
    maybe_write_byte_array_encoding(writer, DataSeries::InsertionBases, encodings.insertion_bases(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::ReferenceSkipLengths, encodings.reference_skip_lengths(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::PaddingLengths, encodings.padding_lengths(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::HardClipLengths, encodings.hard_clip_lengths(), version)?;
    maybe_write_byte_array_encoding(writer, DataSeries::SoftClipBases, encodings.soft_clip_bases(), version)?;
    maybe_write_integer_encoding(writer, DataSeries::MappingQualities, encodings.mapping_qualities(), version)?;
    maybe_write_byte_encoding(writer, DataSeries::Bases, encodings.bases(), version)?;
    maybe_write_byte_encoding(writer, DataSeries::QualityScores, encodings.quality_scores(), version)?;

    Ok(())
}

fn maybe_write_byte_encoding<W>(
    writer: &mut W,
    key: DataSeries,
    encoding: Option<&Encoding<Byte>>,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    if let Some(encoding) = encoding {
        write_key(writer, key)?;
        write_byte_encoding(writer, encoding, version)?;
    }

    Ok(())
}

fn maybe_write_integer_encoding<W>(
    writer: &mut W,
    key: DataSeries,
    encoding: Option<&Encoding<Integer>>,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    if let Some(encoding) = encoding {
        write_key(writer, key)?;
        write_integer_encoding(writer, encoding, version)?;
    }

    Ok(())
}

fn maybe_write_byte_array_encoding<W>(
    writer: &mut W,
    key: DataSeries,
    encoding: Option<&Encoding<ByteArray>>,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    if let Some(encoding) = encoding {
        write_key(writer, key)?;
        write_byte_array_encoding(writer, encoding, version)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_data_series_encodings_len() {
        let encodings = DataSeriesEncodings::default();
        assert_eq!(data_series_encodings_len(&encodings), 0);

        let encodings = DataSeriesEncodings::init(Version::default());
        assert_eq!(data_series_encodings_len(&encodings), 28);

        let encodings = DataSeriesEncodings {
            bam_flags: Some(Encoding::new(Integer::External {
                block_content_id: 1,
            })),
            cram_flags: Some(Encoding::new(Integer::External {
                block_content_id: 2,
            })),
            read_lengths: Some(Encoding::new(Integer::External {
                block_content_id: 4,
            })),
            alignment_starts: Some(Encoding::new(Integer::External {
                block_content_id: 5,
            })),
            read_group_ids: Some(Encoding::new(Integer::External {
                block_content_id: 6,
            })),
            tag_set_ids: Some(Encoding::new(Integer::External {
                block_content_id: 13,
            })),
            ..Default::default()
        };

        assert_eq!(data_series_encodings_len(&encodings), 6);
    }
}
