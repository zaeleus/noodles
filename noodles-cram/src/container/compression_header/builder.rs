use crate::Record;

use super::{
    data_series_encoding_map::{DataSeries, DataSeriesEncodingMap},
    encoding::Encoding,
    preservation_map, tag_encoding_map, CompressionHeader,
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

#[derive(Debug, Default)]
pub struct Builder {
    preservation_map_builder: preservation_map::Builder,
    tag_encoding_map_builder: tag_encoding_map::Builder,
}

impl Builder {
    pub fn update(&mut self, reference_sequence: &[u8], record: &Record) {
        self.preservation_map_builder
            .update(reference_sequence, record);
        self.tag_encoding_map_builder.update(record);
    }

    pub fn build(self) -> CompressionHeader {
        let preservation_map = self.preservation_map_builder.build();

        let mut data_series_encoding_map = DataSeriesEncodingMap::default();

        for (i, &data_series) in DATA_SERIES.iter().enumerate() {
            let block_content_id = (i + 1) as i32;
            // TODO: Select encoding depending on the type of data.
            let encoding = Encoding::External(block_content_id);
            data_series_encoding_map.insert(data_series, encoding);
        }

        data_series_encoding_map.insert(DataSeries::ReadNames, Encoding::ByteArrayStop(0x00, 7));

        let tag_encoding_map = self.tag_encoding_map_builder.build();

        CompressionHeader::new(preservation_map, data_series_encoding_map, tag_encoding_map)
    }
}
