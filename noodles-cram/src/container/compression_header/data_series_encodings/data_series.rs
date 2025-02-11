use std::{error, fmt};

use crate::container::block;

pub static STANDARD_DATA_SERIES: &[DataSeries; 28] = &[
    DataSeries::BamFlags,
    DataSeries::CramFlags,
    DataSeries::ReferenceSequenceIds,
    DataSeries::ReadLengths,
    DataSeries::AlignmentStarts,
    DataSeries::ReadGroupIds,
    DataSeries::Names,
    DataSeries::MateFlags,
    DataSeries::MateReferenceSequenceId,
    DataSeries::MateAlignmentStart,
    DataSeries::TemplateLengths,
    DataSeries::MateDistances,
    DataSeries::TagSetIds,
    DataSeries::FeatureCounts,
    DataSeries::FeatureCodes,
    DataSeries::FeaturePositionDeltas,
    DataSeries::DeletionLengths,
    DataSeries::StretchesOfBases,
    DataSeries::StretchesOfQualityScores,
    DataSeries::BaseSubstitutionCodes,
    DataSeries::InsertionBases,
    DataSeries::ReferenceSkipLengths,
    DataSeries::PaddingLengths,
    DataSeries::HardClipLengths,
    DataSeries::SoftClipBases,
    DataSeries::MappingQualities,
    DataSeries::Bases,
    DataSeries::QualityScores,
];

/// A CRAM container compression header data series.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum DataSeries {
    /// BAM bit flags (`BF`).
    BamFlags,
    /// CRAM bit flags (`CF`).
    CramFlags,
    /// Reference ID (`RI`).
    ReferenceSequenceIds,
    /// Read lengths (`RL`).
    ReadLengths,
    /// In-seq positions (`AP`).
    AlignmentStarts,
    /// Read groups (`RG`).
    ReadGroupIds,
    /// Read names (`RN`).
    Names,
    /// Next mate bit flags (`MF`).
    MateFlags,
    /// Next fragment reference sequence ID (`NS`).
    MateReferenceSequenceId,
    /// Next mate alignment start (`NP`).
    MateAlignmentStart,
    /// Template size (`TS`).
    TemplateLengths,
    /// Distance to next fragment (`NF`).
    MateDistances,
    /// Tag IDs (`TL`).
    TagSetIds,
    /// Number of read features (`FN`).
    FeatureCounts,
    /// Read features codes (`FC`).
    FeatureCodes,
    /// In-read positions (`FP`).
    FeaturePositionDeltas,
    /// Deletion lengths (`DL`).
    DeletionLengths,
    /// Stretches of bases (`BB`).
    StretchesOfBases,
    /// Stretches of quality scores (`QQ`).
    StretchesOfQualityScores,
    /// Base substitution codes (`BS`).
    BaseSubstitutionCodes,
    /// Insertion (`IN`).
    InsertionBases,
    /// Reference skip length (`RS`).
    ReferenceSkipLengths,
    /// Padding (`PD`).
    PaddingLengths,
    /// Hard clip (`HC`).
    HardClipLengths,
    /// Soft clip (`SC`).
    SoftClipBases,
    /// Mapping qualities (`MQ`).
    MappingQualities,
    /// Bases (`BA`).
    Bases,
    /// Quality scores (`QS`).
    QualityScores,
    /// Read tag counts (`TC`).
    ///
    /// This is a legacy CRAM 1.0 data series.
    ReservedTc,
    /// Read tag names and types (`TN`).
    ///
    /// This is a legacy CRAM 1.0 data series.
    ReservedTn,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteArrayError([u8; 2]);

impl error::Error for TryFromByteArrayError {}

impl fmt::Display for TryFromByteArrayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid data series: {:#x?}", self.0)
    }
}

impl TryFrom<[u8; 2]> for DataSeries {
    type Error = TryFromByteArrayError;

    fn try_from(b: [u8; 2]) -> Result<Self, Self::Error> {
        match b {
            [b'B', b'F'] => Ok(Self::BamFlags),
            [b'C', b'F'] => Ok(Self::CramFlags),
            [b'R', b'I'] => Ok(Self::ReferenceSequenceIds),
            [b'R', b'L'] => Ok(Self::ReadLengths),
            [b'A', b'P'] => Ok(Self::AlignmentStarts),
            [b'R', b'G'] => Ok(Self::ReadGroupIds),
            [b'R', b'N'] => Ok(Self::Names),
            [b'M', b'F'] => Ok(Self::MateFlags),
            [b'N', b'S'] => Ok(Self::MateReferenceSequenceId),
            [b'N', b'P'] => Ok(Self::MateAlignmentStart),
            [b'T', b'S'] => Ok(Self::TemplateLengths),
            [b'N', b'F'] => Ok(Self::MateDistances),
            [b'T', b'L'] => Ok(Self::TagSetIds),
            [b'F', b'N'] => Ok(Self::FeatureCounts),
            [b'F', b'C'] => Ok(Self::FeatureCodes),
            [b'F', b'P'] => Ok(Self::FeaturePositionDeltas),
            [b'D', b'L'] => Ok(Self::DeletionLengths),
            [b'B', b'B'] => Ok(Self::StretchesOfBases),
            [b'Q', b'Q'] => Ok(Self::StretchesOfQualityScores),
            [b'B', b'S'] => Ok(Self::BaseSubstitutionCodes),
            [b'I', b'N'] => Ok(Self::InsertionBases),
            [b'R', b'S'] => Ok(Self::ReferenceSkipLengths),
            [b'P', b'D'] => Ok(Self::PaddingLengths),
            [b'H', b'C'] => Ok(Self::HardClipLengths),
            [b'S', b'C'] => Ok(Self::SoftClipBases),
            [b'M', b'Q'] => Ok(Self::MappingQualities),
            [b'B', b'A'] => Ok(Self::Bases),
            [b'Q', b'S'] => Ok(Self::QualityScores),
            [b'T', b'C'] => Ok(Self::ReservedTc),
            [b'T', b'N'] => Ok(Self::ReservedTn),
            _ => Err(TryFromByteArrayError(b)),
        }
    }
}

impl From<DataSeries> for [u8; 2] {
    fn from(data_series: DataSeries) -> Self {
        match data_series {
            DataSeries::BamFlags => [b'B', b'F'],
            DataSeries::CramFlags => [b'C', b'F'],
            DataSeries::ReferenceSequenceIds => [b'R', b'I'],
            DataSeries::ReadLengths => [b'R', b'L'],
            DataSeries::AlignmentStarts => [b'A', b'P'],
            DataSeries::ReadGroupIds => [b'R', b'G'],
            DataSeries::Names => [b'R', b'N'],
            DataSeries::MateFlags => [b'M', b'F'],
            DataSeries::MateReferenceSequenceId => [b'N', b'S'],
            DataSeries::MateAlignmentStart => [b'N', b'P'],
            DataSeries::TemplateLengths => [b'T', b'S'],
            DataSeries::MateDistances => [b'N', b'F'],
            DataSeries::TagSetIds => [b'T', b'L'],
            DataSeries::FeatureCounts => [b'F', b'N'],
            DataSeries::FeatureCodes => [b'F', b'C'],
            DataSeries::FeaturePositionDeltas => [b'F', b'P'],
            DataSeries::DeletionLengths => [b'D', b'L'],
            DataSeries::StretchesOfBases => [b'B', b'B'],
            DataSeries::StretchesOfQualityScores => [b'Q', b'Q'],
            DataSeries::BaseSubstitutionCodes => [b'B', b'S'],
            DataSeries::InsertionBases => [b'I', b'N'],
            DataSeries::ReferenceSkipLengths => [b'R', b'S'],
            DataSeries::PaddingLengths => [b'P', b'D'],
            DataSeries::HardClipLengths => [b'H', b'C'],
            DataSeries::SoftClipBases => [b'S', b'C'],
            DataSeries::MappingQualities => [b'M', b'Q'],
            DataSeries::Bases => [b'B', b'A'],
            DataSeries::QualityScores => [b'Q', b'S'],
            DataSeries::ReservedTc => [b'T', b'C'],
            DataSeries::ReservedTn => [b'T', b'N'],
        }
    }
}

impl From<DataSeries> for block::ContentId {
    fn from(data_series: DataSeries) -> Self {
        match data_series {
            DataSeries::BamFlags => 1,
            DataSeries::CramFlags => 2,
            DataSeries::ReferenceSequenceIds => 3,
            DataSeries::ReadLengths => 4,
            DataSeries::AlignmentStarts => 5,
            DataSeries::ReadGroupIds => 6,
            DataSeries::Names => 7,
            DataSeries::MateFlags => 8,
            DataSeries::MateReferenceSequenceId => 9,
            DataSeries::MateAlignmentStart => 10,
            DataSeries::TemplateLengths => 11,
            DataSeries::MateDistances => 12,
            DataSeries::TagSetIds => 13,
            DataSeries::FeatureCounts => 14,
            DataSeries::FeatureCodes => 15,
            DataSeries::FeaturePositionDeltas => 16,
            DataSeries::DeletionLengths => 17,
            DataSeries::StretchesOfBases => 18,
            DataSeries::StretchesOfQualityScores => 19,
            DataSeries::BaseSubstitutionCodes => 20,
            DataSeries::InsertionBases => 21,
            DataSeries::ReferenceSkipLengths => 22,
            DataSeries::PaddingLengths => 23,
            DataSeries::HardClipLengths => 24,
            DataSeries::SoftClipBases => 25,
            DataSeries::MappingQualities => 26,
            DataSeries::Bases => 27,
            DataSeries::QualityScores => 28,
            DataSeries::ReservedTc => 29,
            DataSeries::ReservedTn => 30,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_byte_array_for_data_series() {
        assert_eq!(DataSeries::try_from([b'B', b'F']), Ok(DataSeries::BamFlags));
        assert_eq!(
            DataSeries::try_from([b'C', b'F']),
            Ok(DataSeries::CramFlags)
        );

        assert_eq!(
            DataSeries::try_from([b'R', b'I']),
            Ok(DataSeries::ReferenceSequenceIds)
        );
        assert_eq!(
            DataSeries::try_from([b'R', b'L']),
            Ok(DataSeries::ReadLengths)
        );
        assert_eq!(
            DataSeries::try_from([b'A', b'P']),
            Ok(DataSeries::AlignmentStarts)
        );
        assert_eq!(
            DataSeries::try_from([b'R', b'G']),
            Ok(DataSeries::ReadGroupIds)
        );
        assert_eq!(DataSeries::try_from([b'R', b'N']), Ok(DataSeries::Names));
        assert_eq!(
            DataSeries::try_from([b'M', b'F']),
            Ok(DataSeries::MateFlags)
        );
        assert_eq!(
            DataSeries::try_from([b'N', b'S']),
            Ok(DataSeries::MateReferenceSequenceId)
        );
        assert_eq!(
            DataSeries::try_from([b'N', b'P']),
            Ok(DataSeries::MateAlignmentStart)
        );
        assert_eq!(
            DataSeries::try_from([b'T', b'S']),
            Ok(DataSeries::TemplateLengths)
        );
        assert_eq!(
            DataSeries::try_from([b'N', b'F']),
            Ok(DataSeries::MateDistances)
        );
        assert_eq!(
            DataSeries::try_from([b'T', b'L']),
            Ok(DataSeries::TagSetIds)
        );
        assert_eq!(
            DataSeries::try_from([b'F', b'N']),
            Ok(DataSeries::FeatureCounts)
        );
        assert_eq!(
            DataSeries::try_from([b'F', b'C']),
            Ok(DataSeries::FeatureCodes)
        );
        assert_eq!(
            DataSeries::try_from([b'F', b'P']),
            Ok(DataSeries::FeaturePositionDeltas)
        );
        assert_eq!(
            DataSeries::try_from([b'D', b'L']),
            Ok(DataSeries::DeletionLengths)
        );
        assert_eq!(
            DataSeries::try_from([b'B', b'B']),
            Ok(DataSeries::StretchesOfBases)
        );
        assert_eq!(
            DataSeries::try_from([b'Q', b'Q']),
            Ok(DataSeries::StretchesOfQualityScores)
        );
        assert_eq!(
            DataSeries::try_from([b'B', b'S']),
            Ok(DataSeries::BaseSubstitutionCodes)
        );
        assert_eq!(
            DataSeries::try_from([b'I', b'N']),
            Ok(DataSeries::InsertionBases)
        );
        assert_eq!(
            DataSeries::try_from([b'R', b'S']),
            Ok(DataSeries::ReferenceSkipLengths)
        );
        assert_eq!(
            DataSeries::try_from([b'P', b'D']),
            Ok(DataSeries::PaddingLengths)
        );
        assert_eq!(
            DataSeries::try_from([b'H', b'C']),
            Ok(DataSeries::HardClipLengths)
        );
        assert_eq!(
            DataSeries::try_from([b'S', b'C']),
            Ok(DataSeries::SoftClipBases)
        );
        assert_eq!(
            DataSeries::try_from([b'M', b'Q']),
            Ok(DataSeries::MappingQualities)
        );
        assert_eq!(DataSeries::try_from([b'B', b'A']), Ok(DataSeries::Bases));
        assert_eq!(
            DataSeries::try_from([b'Q', b'S']),
            Ok(DataSeries::QualityScores)
        );
        assert_eq!(
            DataSeries::try_from([b'T', b'N']),
            Ok(DataSeries::ReservedTn)
        );
        assert_eq!(
            DataSeries::try_from([b'T', b'C']),
            Ok(DataSeries::ReservedTc)
        );

        assert_eq!(
            DataSeries::try_from([b'X', b'Y']),
            Err(TryFromByteArrayError([b'X', b'Y']))
        );
    }

    #[test]
    fn test_from_data_series_for_u8_array() {
        assert_eq!(<[u8; 2]>::from(DataSeries::BamFlags), [b'B', b'F']);
        assert_eq!(<[u8; 2]>::from(DataSeries::CramFlags), [b'C', b'F']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::ReferenceSequenceIds),
            [b'R', b'I']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::ReadLengths), [b'R', b'L']);
        assert_eq!(<[u8; 2]>::from(DataSeries::AlignmentStarts), [b'A', b'P']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReadGroupIds), [b'R', b'G']);
        assert_eq!(<[u8; 2]>::from(DataSeries::Names), [b'R', b'N']);
        assert_eq!(<[u8; 2]>::from(DataSeries::MateFlags), [b'M', b'F']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::MateReferenceSequenceId),
            [b'N', b'S']
        );
        assert_eq!(
            <[u8; 2]>::from(DataSeries::MateAlignmentStart),
            [b'N', b'P']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::TemplateLengths), [b'T', b'S']);
        assert_eq!(<[u8; 2]>::from(DataSeries::MateDistances), [b'N', b'F']);
        assert_eq!(<[u8; 2]>::from(DataSeries::TagSetIds), [b'T', b'L']);
        assert_eq!(<[u8; 2]>::from(DataSeries::FeatureCounts), [b'F', b'N']);
        assert_eq!(<[u8; 2]>::from(DataSeries::FeatureCodes), [b'F', b'C']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::FeaturePositionDeltas),
            [b'F', b'P']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::DeletionLengths), [b'D', b'L']);
        assert_eq!(<[u8; 2]>::from(DataSeries::StretchesOfBases), [b'B', b'B']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::StretchesOfQualityScores),
            [b'Q', b'Q']
        );
        assert_eq!(
            <[u8; 2]>::from(DataSeries::BaseSubstitutionCodes),
            [b'B', b'S']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::InsertionBases), [b'I', b'N']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::ReferenceSkipLengths),
            [b'R', b'S']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::PaddingLengths), [b'P', b'D']);
        assert_eq!(<[u8; 2]>::from(DataSeries::HardClipLengths), [b'H', b'C']);
        assert_eq!(<[u8; 2]>::from(DataSeries::SoftClipBases), [b'S', b'C']);
        assert_eq!(<[u8; 2]>::from(DataSeries::MappingQualities), [b'M', b'Q']);
        assert_eq!(<[u8; 2]>::from(DataSeries::Bases), [b'B', b'A']);
        assert_eq!(<[u8; 2]>::from(DataSeries::QualityScores), [b'Q', b'S']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReservedTn), [b'T', b'N']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReservedTc), [b'T', b'C']);
    }
}
