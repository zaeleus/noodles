use std::io::{self, Write};

use crate::{
    container::compression_header::{data_series_encoding_map::DataSeries, DataSeriesEncodingMap},
    num::{write_itf8, Itf8},
    writer::encoding::write_encoding,
};

const DATA_SERIES: [DataSeries; DataSeries::LEN] = [
    DataSeries::BamBitFlags,
    DataSeries::CramBitFlags,
    DataSeries::ReferenceId,
    DataSeries::ReadLengths,
    DataSeries::InSeqPositions,
    DataSeries::ReadGroups,
    DataSeries::ReadNames,
    DataSeries::NextMateBitFlags,
    DataSeries::NextFragmentReferenceSequenceId,
    DataSeries::NextMateAlignmentStart,
    DataSeries::TemplateSize,
    DataSeries::DistanceToNextFragment,
    DataSeries::TagIds,
    DataSeries::NumberOfReadFeatures,
    DataSeries::ReadFeaturesCodes,
    DataSeries::InReadPositions,
    DataSeries::DeletionLengths,
    DataSeries::StretchesOfBases,
    DataSeries::StretchesOfQualityScores,
    DataSeries::BaseSubstitutionCodes,
    DataSeries::Insertion,
    DataSeries::ReferenceSkipLength,
    DataSeries::Padding,
    DataSeries::HardClip,
    DataSeries::SoftClip,
    DataSeries::MappingQualities,
    DataSeries::Bases,
    DataSeries::QualityScores,
];

#[allow(dead_code)]
pub fn write_data_series_encoding_map<W>(
    writer: &mut W,
    data_series_encoding_map: &DataSeriesEncodingMap,
) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    let mut map_len = 0;

    for data_series in &DATA_SERIES {
        if let Some(encoding) = data_series_encoding_map.get(data_series) {
            write_key(&mut buf, *data_series)?;
            write_encoding(&mut buf, encoding)?;
            map_len += 1;
        }
    }

    // FIXME: usize => Itf8 cast
    let data_len = buf.len() as Itf8;
    write_itf8(writer, data_len)?;

    write_itf8(writer, map_len)?;

    writer.write_all(&buf)
}

fn write_key<W>(writer: &mut W, key: DataSeries) -> io::Result<()>
where
    W: Write,
{
    let data = <[u8; 2]>::from(key);
    writer.write_all(&data)
}
