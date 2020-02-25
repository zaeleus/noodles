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
    BamBitFlags,
    CramBitFlags,
    ReferenceId,
    ReadLengths,
    InSeqPositions,
    ReadGroups,
    ReadNames,
    NextMateBitFlags,
    NextFragmentReferenceSequenceId,
    NextMateAlignmentStart,
    TemplateSize,
    DistanceToNextFragment,
    TagIds,
    NumberOfReadFeatures,
    ReadFeaturesCodes,
    InReadPositions,
    DeletionLengths,
    StretchesOfBases,
    StretchesOfQualityScores,
    BaseSubstitutionCodes,
    Insertion,
    ReferenceSkipLength,
    Padding,
    HardClip,
    SoftClip,
    MappingQualities,
    Bases,
    QualityScores,
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

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use super::*;

    #[test]
    fn test_try_from() {
        assert_eq!(
            DataSeries::try_from(&b"BF"[..]),
            Ok(DataSeries::BamBitFlags)
        );
        assert_eq!(
            DataSeries::try_from(&b"CF"[..]),
            Ok(DataSeries::CramBitFlags)
        );
        assert_eq!(
            DataSeries::try_from(&b"XY"[..]),
            Err(TryFromByteSliceError(vec![0x58, 0x59]))
        );
    }
}
