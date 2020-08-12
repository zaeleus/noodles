pub mod data_series_encoding_map;
pub mod encoding;
pub mod preservation_map;
mod tag_encoding_map;

pub use self::{
    data_series_encoding_map::DataSeriesEncodingMap,
    encoding::Encoding,
    preservation_map::{PreservationMap, SubstitutionMatrix, TagIdsDictionary},
    tag_encoding_map::TagEncodingMap,
};

use std::{convert::TryFrom, io};

use crate::{reader::compression_header::read_compression_header, Record};

use super::Block;

#[derive(Debug)]
pub struct CompressionHeader {
    preservation_map: PreservationMap,
    data_series_encoding_map: DataSeriesEncodingMap,
    tag_encoding_map: TagEncodingMap,
}

impl CompressionHeader {
    pub fn from_records(reference_sequence: &[u8], records: &[Record]) -> Self {
        use data_series_encoding_map::DataSeries;

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

        let mut preservation_map_builder = PreservationMap::builder(reference_sequence);

        for record in records {
            preservation_map_builder.update(record);
        }

        let preservation_map = preservation_map_builder.build();

        let mut data_series_encoding_map = DataSeriesEncodingMap::default();

        for (i, &data_series) in DATA_SERIES.iter().enumerate() {
            let block_content_id = (i + 1) as i32;
            // TODO: Select encoding depending on the type of data.
            let encoding = Encoding::External(block_content_id);
            data_series_encoding_map.insert(data_series, encoding);
        }

        data_series_encoding_map.insert(DataSeries::ReadNames, Encoding::ByteArrayStop(0x00, 7));

        let tag_encoding_map = TagEncodingMap::from_records(records);

        Self::new(preservation_map, data_series_encoding_map, tag_encoding_map)
    }

    pub fn new(
        preservation_map: PreservationMap,
        data_series_encoding_map: DataSeriesEncodingMap,
        tag_encoding_map: TagEncodingMap,
    ) -> Self {
        Self {
            preservation_map,
            data_series_encoding_map,
            tag_encoding_map,
        }
    }

    pub fn preservation_map(&self) -> &PreservationMap {
        &self.preservation_map
    }

    pub fn data_series_encoding_map(&self) -> &DataSeriesEncodingMap {
        &self.data_series_encoding_map
    }

    pub fn tag_encoding_map(&self) -> &TagEncodingMap {
        &self.tag_encoding_map
    }
}

impl TryFrom<&Block> for CompressionHeader {
    type Error = io::Error;

    fn try_from(block: &Block) -> Result<Self, Self::Error> {
        let data = block.decompressed_data();
        let mut reader = &data[..];
        read_compression_header(&mut reader)
    }
}
