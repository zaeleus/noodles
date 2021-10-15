use std::{convert::TryFrom, error, fmt};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum DataSeries {
    /// BAM bit flags (`BF`).
    BamBitFlags,
    /// CRAM bit flags (`CF`).
    CramBitFlags,
    /// Reference ID (`RI`).
    ReferenceId,
    /// Read lengths (`RL`).
    ReadLengths,
    /// In-seq positions (`AP`).
    InSeqPositions,
    /// Read groups (`RG`).
    ReadGroups,
    /// Read names (`RN`).
    ReadNames,
    /// Next mate bit flags (`MF`).
    NextMateBitFlags,
    /// Next fragment reference sequence ID (`NS`).
    NextFragmentReferenceSequenceId,
    /// Next mate alignment start (`NP`).
    NextMateAlignmentStart,
    /// Template size (`TS`).
    TemplateSize,
    /// Distance to next fragment (`NF`).
    DistanceToNextFragment,
    /// Tag IDs (`TL`).
    TagIds,
    /// Number of read features (`FN`).
    NumberOfReadFeatures,
    /// Read features codes (`FC`).
    ReadFeaturesCodes,
    /// In-read positions (`FP`).
    InReadPositions,
    /// Deletion lengths (`DL`).
    DeletionLengths,
    /// Stretches of bases (`BB`).
    StretchesOfBases,
    /// Stretches of quality scores (`QQ`).
    StretchesOfQualityScores,
    /// Base substitution codes (`BS`).
    BaseSubstitutionCodes,
    /// Insertion (`IN`).
    Insertion,
    /// Reference skip length (`RS`).
    ReferenceSkipLength,
    /// Pading (`PD`).
    Padding,
    /// Hard clip (`HC`).
    HardClip,
    /// Soft clip (`SC`).
    SoftClip,
    /// Mapping qualities (`MQ`).
    MappingQualities,
    /// Bases (`BA`).
    Bases,
    /// Quality scores (`QS`).
    QualityScores,
    /// (`TC`).
    ///
    /// Legacy CRAM 1.0 data series.
    ReservedTc,
    /// (`TN`).
    ///
    /// Legacy CRAM 1.0 data series.
    ReservedTn,
}

impl DataSeries {
    /// The number of data series variants.
    pub(crate) const LEN: usize = 28;
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
            [b'B', b'F'] => Ok(Self::BamBitFlags),
            [b'C', b'F'] => Ok(Self::CramBitFlags),
            [b'R', b'I'] => Ok(Self::ReferenceId),
            [b'R', b'L'] => Ok(Self::ReadLengths),
            [b'A', b'P'] => Ok(Self::InSeqPositions),
            [b'R', b'G'] => Ok(Self::ReadGroups),
            [b'R', b'N'] => Ok(Self::ReadNames),
            [b'M', b'F'] => Ok(Self::NextMateBitFlags),
            [b'N', b'S'] => Ok(Self::NextFragmentReferenceSequenceId),
            [b'N', b'P'] => Ok(Self::NextMateAlignmentStart),
            [b'T', b'S'] => Ok(Self::TemplateSize),
            [b'N', b'F'] => Ok(Self::DistanceToNextFragment),
            [b'T', b'L'] => Ok(Self::TagIds),
            [b'F', b'N'] => Ok(Self::NumberOfReadFeatures),
            [b'F', b'C'] => Ok(Self::ReadFeaturesCodes),
            [b'F', b'P'] => Ok(Self::InReadPositions),
            [b'D', b'L'] => Ok(Self::DeletionLengths),
            [b'B', b'B'] => Ok(Self::StretchesOfBases),
            [b'Q', b'Q'] => Ok(Self::StretchesOfQualityScores),
            [b'B', b'S'] => Ok(Self::BaseSubstitutionCodes),
            [b'I', b'N'] => Ok(Self::Insertion),
            [b'R', b'S'] => Ok(Self::ReferenceSkipLength),
            [b'P', b'D'] => Ok(Self::Padding),
            [b'H', b'C'] => Ok(Self::HardClip),
            [b'S', b'C'] => Ok(Self::SoftClip),
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
            DataSeries::BamBitFlags => [b'B', b'F'],
            DataSeries::CramBitFlags => [b'C', b'F'],
            DataSeries::ReferenceId => [b'R', b'I'],
            DataSeries::ReadLengths => [b'R', b'L'],
            DataSeries::InSeqPositions => [b'A', b'P'],
            DataSeries::ReadGroups => [b'R', b'G'],
            DataSeries::ReadNames => [b'R', b'N'],
            DataSeries::NextMateBitFlags => [b'M', b'F'],
            DataSeries::NextFragmentReferenceSequenceId => [b'N', b'S'],
            DataSeries::NextMateAlignmentStart => [b'N', b'P'],
            DataSeries::TemplateSize => [b'T', b'S'],
            DataSeries::DistanceToNextFragment => [b'N', b'F'],
            DataSeries::TagIds => [b'T', b'L'],
            DataSeries::NumberOfReadFeatures => [b'F', b'N'],
            DataSeries::ReadFeaturesCodes => [b'F', b'C'],
            DataSeries::InReadPositions => [b'F', b'P'],
            DataSeries::DeletionLengths => [b'D', b'L'],
            DataSeries::StretchesOfBases => [b'B', b'B'],
            DataSeries::StretchesOfQualityScores => [b'Q', b'Q'],
            DataSeries::BaseSubstitutionCodes => [b'B', b'S'],
            DataSeries::Insertion => [b'I', b'N'],
            DataSeries::ReferenceSkipLength => [b'R', b'S'],
            DataSeries::Padding => [b'P', b'D'],
            DataSeries::HardClip => [b'H', b'C'],
            DataSeries::SoftClip => [b'S', b'C'],
            DataSeries::MappingQualities => [b'M', b'Q'],
            DataSeries::Bases => [b'B', b'A'],
            DataSeries::QualityScores => [b'Q', b'S'],
            DataSeries::ReservedTc => [b'T', b'C'],
            DataSeries::ReservedTn => [b'T', b'N'],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_byte_array_for_data_series() {
        assert_eq!(
            DataSeries::try_from([b'B', b'F']),
            Ok(DataSeries::BamBitFlags)
        );
        assert_eq!(
            DataSeries::try_from([b'C', b'F']),
            Ok(DataSeries::CramBitFlags)
        );

        assert_eq!(
            DataSeries::try_from([b'R', b'I']),
            Ok(DataSeries::ReferenceId)
        );
        assert_eq!(
            DataSeries::try_from([b'R', b'L']),
            Ok(DataSeries::ReadLengths)
        );
        assert_eq!(
            DataSeries::try_from([b'A', b'P']),
            Ok(DataSeries::InSeqPositions)
        );
        assert_eq!(
            DataSeries::try_from([b'R', b'G']),
            Ok(DataSeries::ReadGroups)
        );
        assert_eq!(
            DataSeries::try_from([b'R', b'N']),
            Ok(DataSeries::ReadNames)
        );
        assert_eq!(
            DataSeries::try_from([b'M', b'F']),
            Ok(DataSeries::NextMateBitFlags)
        );
        assert_eq!(
            DataSeries::try_from([b'N', b'S']),
            Ok(DataSeries::NextFragmentReferenceSequenceId)
        );
        assert_eq!(
            DataSeries::try_from([b'N', b'P']),
            Ok(DataSeries::NextMateAlignmentStart)
        );
        assert_eq!(
            DataSeries::try_from([b'T', b'S']),
            Ok(DataSeries::TemplateSize)
        );
        assert_eq!(
            DataSeries::try_from([b'N', b'F']),
            Ok(DataSeries::DistanceToNextFragment)
        );
        assert_eq!(DataSeries::try_from([b'T', b'L']), Ok(DataSeries::TagIds));
        assert_eq!(
            DataSeries::try_from([b'F', b'N']),
            Ok(DataSeries::NumberOfReadFeatures)
        );
        assert_eq!(
            DataSeries::try_from([b'F', b'C']),
            Ok(DataSeries::ReadFeaturesCodes)
        );
        assert_eq!(
            DataSeries::try_from([b'F', b'P']),
            Ok(DataSeries::InReadPositions)
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
            Ok(DataSeries::Insertion)
        );
        assert_eq!(
            DataSeries::try_from([b'R', b'S']),
            Ok(DataSeries::ReferenceSkipLength)
        );
        assert_eq!(DataSeries::try_from([b'P', b'D']), Ok(DataSeries::Padding));
        assert_eq!(DataSeries::try_from([b'H', b'C']), Ok(DataSeries::HardClip));
        assert_eq!(DataSeries::try_from([b'S', b'C']), Ok(DataSeries::SoftClip));
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
        assert_eq!(<[u8; 2]>::from(DataSeries::BamBitFlags), [b'B', b'F']);
        assert_eq!(<[u8; 2]>::from(DataSeries::CramBitFlags), [b'C', b'F']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReferenceId), [b'R', b'I']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReadLengths), [b'R', b'L']);
        assert_eq!(<[u8; 2]>::from(DataSeries::InSeqPositions), [b'A', b'P']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReadGroups), [b'R', b'G']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReadNames), [b'R', b'N']);
        assert_eq!(<[u8; 2]>::from(DataSeries::NextMateBitFlags), [b'M', b'F']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::NextFragmentReferenceSequenceId),
            [b'N', b'S']
        );
        assert_eq!(
            <[u8; 2]>::from(DataSeries::NextMateAlignmentStart),
            [b'N', b'P']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::TemplateSize), [b'T', b'S']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::DistanceToNextFragment),
            [b'N', b'F']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::TagIds), [b'T', b'L']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::NumberOfReadFeatures),
            [b'F', b'N']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::ReadFeaturesCodes), [b'F', b'C']);
        assert_eq!(<[u8; 2]>::from(DataSeries::InReadPositions), [b'F', b'P']);
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
        assert_eq!(<[u8; 2]>::from(DataSeries::Insertion), [b'I', b'N']);
        assert_eq!(
            <[u8; 2]>::from(DataSeries::ReferenceSkipLength),
            [b'R', b'S']
        );
        assert_eq!(<[u8; 2]>::from(DataSeries::Padding), [b'P', b'D']);
        assert_eq!(<[u8; 2]>::from(DataSeries::HardClip), [b'H', b'C']);
        assert_eq!(<[u8; 2]>::from(DataSeries::SoftClip), [b'S', b'C']);
        assert_eq!(<[u8; 2]>::from(DataSeries::MappingQualities), [b'M', b'Q']);
        assert_eq!(<[u8; 2]>::from(DataSeries::Bases), [b'B', b'A']);
        assert_eq!(<[u8; 2]>::from(DataSeries::QualityScores), [b'Q', b'S']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReservedTn), [b'T', b'N']);
        assert_eq!(<[u8; 2]>::from(DataSeries::ReservedTc), [b'T', b'C']);
    }
}
