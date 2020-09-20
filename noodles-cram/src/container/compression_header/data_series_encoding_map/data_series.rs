use std::{convert::TryFrom, error, fmt};

#[derive(Debug, Eq, PartialEq)]
pub struct TryFromByteSliceError(Vec<u8>);

impl fmt::Display for TryFromByteSliceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid data series: {:#x?}", self.0)
    }
}

impl error::Error for TryFromByteSliceError {}

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
}

impl DataSeries {
    /// The number of data series variants.
    pub(crate) const LEN: usize = 28;
}

impl TryFrom<&[u8]> for DataSeries {
    type Error = TryFromByteSliceError;

    fn try_from(b: &[u8]) -> Result<Self, Self::Error> {
        match b {
            b"BF" => Ok(Self::BamBitFlags),
            b"CF" => Ok(Self::CramBitFlags),
            b"RI" => Ok(Self::ReferenceId),
            b"RL" => Ok(Self::ReadLengths),
            b"AP" => Ok(Self::InSeqPositions),
            b"RG" => Ok(Self::ReadGroups),
            b"RN" => Ok(Self::ReadNames),
            b"MF" => Ok(Self::NextMateBitFlags),
            b"NS" => Ok(Self::NextFragmentReferenceSequenceId),
            b"NP" => Ok(Self::NextMateAlignmentStart),
            b"TS" => Ok(Self::TemplateSize),
            b"NF" => Ok(Self::DistanceToNextFragment),
            b"TL" => Ok(Self::TagIds),
            b"FN" => Ok(Self::NumberOfReadFeatures),
            b"FC" => Ok(Self::ReadFeaturesCodes),
            b"FP" => Ok(Self::InReadPositions),
            b"DL" => Ok(Self::DeletionLengths),
            b"BB" => Ok(Self::StretchesOfBases),
            b"QQ" => Ok(Self::StretchesOfQualityScores),
            b"BS" => Ok(Self::BaseSubstitutionCodes),
            b"IN" => Ok(Self::Insertion),
            b"RS" => Ok(Self::ReferenceSkipLength),
            b"PD" => Ok(Self::Padding),
            b"HC" => Ok(Self::HardClip),
            b"SC" => Ok(Self::SoftClip),
            b"MQ" => Ok(Self::MappingQualities),
            b"BA" => Ok(Self::Bases),
            b"QS" => Ok(Self::QualityScores),
            _ => Err(TryFromByteSliceError(b.to_vec())),
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
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_byte_slice_for_data_series() {
        assert_eq!(
            DataSeries::try_from(&b"BF"[..]),
            Ok(DataSeries::BamBitFlags)
        );
        assert_eq!(
            DataSeries::try_from(&b"CF"[..]),
            Ok(DataSeries::CramBitFlags)
        );

        assert_eq!(
            DataSeries::try_from(&b"RI"[..]),
            Ok(DataSeries::ReferenceId)
        );
        assert_eq!(
            DataSeries::try_from(&b"RL"[..]),
            Ok(DataSeries::ReadLengths)
        );
        assert_eq!(
            DataSeries::try_from(&b"AP"[..]),
            Ok(DataSeries::InSeqPositions)
        );
        assert_eq!(DataSeries::try_from(&b"RG"[..]), Ok(DataSeries::ReadGroups));
        assert_eq!(DataSeries::try_from(&b"RN"[..]), Ok(DataSeries::ReadNames));
        assert_eq!(
            DataSeries::try_from(&b"MF"[..]),
            Ok(DataSeries::NextMateBitFlags)
        );
        assert_eq!(
            DataSeries::try_from(&b"NS"[..]),
            Ok(DataSeries::NextFragmentReferenceSequenceId)
        );
        assert_eq!(
            DataSeries::try_from(&b"NP"[..]),
            Ok(DataSeries::NextMateAlignmentStart)
        );
        assert_eq!(
            DataSeries::try_from(&b"TS"[..]),
            Ok(DataSeries::TemplateSize)
        );
        assert_eq!(
            DataSeries::try_from(&b"NF"[..]),
            Ok(DataSeries::DistanceToNextFragment)
        );
        assert_eq!(DataSeries::try_from(&b"TL"[..]), Ok(DataSeries::TagIds));
        assert_eq!(
            DataSeries::try_from(&b"FN"[..]),
            Ok(DataSeries::NumberOfReadFeatures)
        );
        assert_eq!(
            DataSeries::try_from(&b"FC"[..]),
            Ok(DataSeries::ReadFeaturesCodes)
        );
        assert_eq!(
            DataSeries::try_from(&b"FP"[..]),
            Ok(DataSeries::InReadPositions)
        );
        assert_eq!(
            DataSeries::try_from(&b"DL"[..]),
            Ok(DataSeries::DeletionLengths)
        );
        assert_eq!(
            DataSeries::try_from(&b"BB"[..]),
            Ok(DataSeries::StretchesOfBases)
        );
        assert_eq!(
            DataSeries::try_from(&b"QQ"[..]),
            Ok(DataSeries::StretchesOfQualityScores)
        );
        assert_eq!(
            DataSeries::try_from(&b"BS"[..]),
            Ok(DataSeries::BaseSubstitutionCodes)
        );
        assert_eq!(DataSeries::try_from(&b"IN"[..]), Ok(DataSeries::Insertion));
        assert_eq!(
            DataSeries::try_from(&b"RS"[..]),
            Ok(DataSeries::ReferenceSkipLength)
        );
        assert_eq!(DataSeries::try_from(&b"PD"[..]), Ok(DataSeries::Padding));
        assert_eq!(DataSeries::try_from(&b"HC"[..]), Ok(DataSeries::HardClip));
        assert_eq!(DataSeries::try_from(&b"SC"[..]), Ok(DataSeries::SoftClip));
        assert_eq!(
            DataSeries::try_from(&b"MQ"[..]),
            Ok(DataSeries::MappingQualities)
        );
        assert_eq!(DataSeries::try_from(&b"BA"[..]), Ok(DataSeries::Bases));
        assert_eq!(
            DataSeries::try_from(&b"QS"[..]),
            Ok(DataSeries::QualityScores)
        );

        assert_eq!(
            DataSeries::try_from(&b"XY"[..]),
            Err(TryFromByteSliceError(vec![0x58, 0x59]))
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
    }
}
