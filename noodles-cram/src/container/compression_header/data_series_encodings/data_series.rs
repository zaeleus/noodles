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
    DataSeries::MateReferenceSequenceIds,
    DataSeries::MateAlignmentStarts,
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
    MateReferenceSequenceIds,
    /// Next mate alignment start (`NP`).
    MateAlignmentStarts,
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
            [b'N', b'S'] => Ok(Self::MateReferenceSequenceIds),
            [b'N', b'P'] => Ok(Self::MateAlignmentStarts),
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
            DataSeries::BamFlags => *b"BF",
            DataSeries::CramFlags => *b"CF",
            DataSeries::ReferenceSequenceIds => *b"RI",
            DataSeries::ReadLengths => *b"RL",
            DataSeries::AlignmentStarts => *b"AP",
            DataSeries::ReadGroupIds => *b"RG",
            DataSeries::Names => *b"RN",
            DataSeries::MateFlags => *b"MF",
            DataSeries::MateReferenceSequenceIds => *b"NS",
            DataSeries::MateAlignmentStarts => *b"NP",
            DataSeries::TemplateLengths => *b"TS",
            DataSeries::MateDistances => *b"NF",
            DataSeries::TagSetIds => *b"TL",
            DataSeries::FeatureCounts => *b"FN",
            DataSeries::FeatureCodes => *b"FC",
            DataSeries::FeaturePositionDeltas => *b"FP",
            DataSeries::DeletionLengths => *b"DL",
            DataSeries::StretchesOfBases => *b"BB",
            DataSeries::StretchesOfQualityScores => *b"QQ",
            DataSeries::BaseSubstitutionCodes => *b"BS",
            DataSeries::InsertionBases => *b"IN",
            DataSeries::ReferenceSkipLengths => *b"RS",
            DataSeries::PaddingLengths => *b"PD",
            DataSeries::HardClipLengths => *b"HC",
            DataSeries::SoftClipBases => *b"SC",
            DataSeries::MappingQualities => *b"MQ",
            DataSeries::Bases => *b"BA",
            DataSeries::QualityScores => *b"QS",
            DataSeries::ReservedTc => *b"TC",
            DataSeries::ReservedTn => *b"TN",
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
            DataSeries::MateReferenceSequenceIds => 9,
            DataSeries::MateAlignmentStarts => 10,
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
        assert_eq!(DataSeries::try_from(*b"BF"), Ok(DataSeries::BamFlags));
        assert_eq!(DataSeries::try_from(*b"CF"), Ok(DataSeries::CramFlags));

        assert_eq!(
            DataSeries::try_from(*b"RI"),
            Ok(DataSeries::ReferenceSequenceIds)
        );
        assert_eq!(DataSeries::try_from(*b"RL"), Ok(DataSeries::ReadLengths));
        assert_eq!(
            DataSeries::try_from(*b"AP"),
            Ok(DataSeries::AlignmentStarts)
        );
        assert_eq!(DataSeries::try_from(*b"RG"), Ok(DataSeries::ReadGroupIds));
        assert_eq!(DataSeries::try_from(*b"RN"), Ok(DataSeries::Names));
        assert_eq!(DataSeries::try_from(*b"MF"), Ok(DataSeries::MateFlags));
        assert_eq!(
            DataSeries::try_from(*b"NS"),
            Ok(DataSeries::MateReferenceSequenceIds)
        );
        assert_eq!(
            DataSeries::try_from(*b"NP"),
            Ok(DataSeries::MateAlignmentStarts)
        );
        assert_eq!(
            DataSeries::try_from(*b"TS"),
            Ok(DataSeries::TemplateLengths)
        );
        assert_eq!(DataSeries::try_from(*b"NF"), Ok(DataSeries::MateDistances));
        assert_eq!(DataSeries::try_from(*b"TL"), Ok(DataSeries::TagSetIds));
        assert_eq!(DataSeries::try_from(*b"FN"), Ok(DataSeries::FeatureCounts));
        assert_eq!(DataSeries::try_from(*b"FC"), Ok(DataSeries::FeatureCodes));
        assert_eq!(
            DataSeries::try_from(*b"FP"),
            Ok(DataSeries::FeaturePositionDeltas)
        );
        assert_eq!(
            DataSeries::try_from(*b"DL"),
            Ok(DataSeries::DeletionLengths)
        );
        assert_eq!(
            DataSeries::try_from(*b"BB"),
            Ok(DataSeries::StretchesOfBases)
        );
        assert_eq!(
            DataSeries::try_from(*b"QQ"),
            Ok(DataSeries::StretchesOfQualityScores)
        );
        assert_eq!(
            DataSeries::try_from(*b"BS"),
            Ok(DataSeries::BaseSubstitutionCodes)
        );
        assert_eq!(DataSeries::try_from(*b"IN"), Ok(DataSeries::InsertionBases));
        assert_eq!(
            DataSeries::try_from(*b"RS"),
            Ok(DataSeries::ReferenceSkipLengths)
        );
        assert_eq!(DataSeries::try_from(*b"PD"), Ok(DataSeries::PaddingLengths));
        assert_eq!(
            DataSeries::try_from(*b"HC"),
            Ok(DataSeries::HardClipLengths)
        );
        assert_eq!(DataSeries::try_from(*b"SC"), Ok(DataSeries::SoftClipBases));
        assert_eq!(
            DataSeries::try_from(*b"MQ"),
            Ok(DataSeries::MappingQualities)
        );
        assert_eq!(DataSeries::try_from(*b"BA"), Ok(DataSeries::Bases));
        assert_eq!(DataSeries::try_from(*b"QS"), Ok(DataSeries::QualityScores));
        assert_eq!(DataSeries::try_from(*b"TN"), Ok(DataSeries::ReservedTn));
        assert_eq!(DataSeries::try_from(*b"TC"), Ok(DataSeries::ReservedTc));

        assert_eq!(
            DataSeries::try_from(*b"XY"),
            Err(TryFromByteArrayError(*b"XY"))
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
            <[u8; 2]>::from(DataSeries::MateReferenceSequenceIds),
            [b'N', b'S']
        );
        assert_eq!(
            <[u8; 2]>::from(DataSeries::MateAlignmentStarts),
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
