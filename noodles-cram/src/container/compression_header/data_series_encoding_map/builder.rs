use crate::container::compression_header::Encoding;

use super::{DataSeries, DataSeriesEncodingMap};

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

pub struct Builder {
    map: [Option<Encoding>; DataSeries::LEN],
}

impl Builder {
    pub fn build(self) -> DataSeriesEncodingMap {
        DataSeriesEncodingMap { map: self.map }
    }
}

impl Default for Builder {
    fn default() -> Self {
        let mut map: [Option<Encoding>; DataSeries::LEN] = Default::default();

        for (i, &data_series) in DATA_SERIES.iter().enumerate() {
            let block_content_id = (i + 1) as i32;
            let encoding = Encoding::External(block_content_id);
            map[data_series as usize] = Some(encoding);
        }

        map[DataSeries::ReadNames as usize] = Some(Encoding::ByteArrayStop(0x00, 7));
        map[DataSeries::StretchesOfBases as usize] = Some(Encoding::ByteArrayStop(0x00, 18));

        map[DataSeries::StretchesOfQualityScores as usize] = Some(Encoding::ByteArrayLen(
            Box::new(Encoding::External(19)),
            Box::new(Encoding::External(19)),
        ));

        map[DataSeries::Insertion as usize] = Some(Encoding::ByteArrayStop(0x00, 21));
        map[DataSeries::SoftClip as usize] = Some(Encoding::ByteArrayStop(0x00, 25));

        Self { map }
    }
}
