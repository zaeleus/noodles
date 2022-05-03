//! SAM record data field tag.

use std::{error, fmt, fmt::Write, str::FromStr};

const LENGTH: usize = 2;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[doc(hidden)]
pub struct Other([u8; LENGTH]);

/// A SAM record data field tag.
///
/// Standard tags are defined in "Sequence Alignment/Map Optional Fields Specification"
/// (2020-05-29).
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Tag {
    /// (`AM`).
    MinMappingQuality,
    /// (`AS`).
    AlignmentScore,
    /// (`BC`).
    SampleBarcodeSequence,
    /// (`BQ`).
    BaseAlignmentQualityOffsets,
    /// (`BZ`).
    OriginalUmiQualityScores,
    /// (`CB`).
    CellBarcodeId,
    /// (`CC`).
    NextHitReferenceSequenceName,
    /// (`CG`).
    Cigar,
    /// (`CM`).
    ColorEditDistance,
    /// (`CO`).
    Comment,
    /// (`CP`).
    NextHitPosition,
    /// (`CQ`).
    ColarQualityScores,
    /// (`CR`).
    CellBarcodeSequence,
    /// (`CS`).
    ColorSequence,
    /// (`CT`).
    CompleteReadAnnotations,
    /// (`CY`).
    CellBarcodeQualityScores,
    /// (`E2`).
    NextHitSequence,
    /// (`FI`).
    SegmentIndex,
    /// (`FS`).
    SegmentSuffix,
    /// (`FZ`).
    AlternativeSequence,
    /// (`GC`).
    ReservedGc,
    /// (`GQ`).
    ReservedGq,
    /// (`GS`).
    ReservedGs,
    /// (`H0`).
    PerfectHitCount,
    /// (`H1`).
    OneDifferenceHitCount,
    /// (`H2`).
    TwoDifferenceHitCount,
    /// (`HI`).
    HitIndex,
    /// (`IH`).
    TotalHitCount,
    /// (`LB`).
    Library,
    /// (`MC`).
    MateCigar,
    /// (`MD`).
    MismatchedPositions,
    /// (`MF`).
    ReservedMf,
    /// (`MI`).
    UmiId,
    /// (`ML`).
    BaseModificationProbabilities,
    /// (`MM`).
    BaseModifications,
    /// (`MQ`).
    MateMappingQuality,
    /// (`NH`).
    AlignmentHitCount,
    /// (`NM`).
    EditDistance,
    /// (`OA`).
    OriginalAlignment,
    /// (`OC`).
    OriginalCigar,
    /// (`OP`).
    OriginalPosition,
    /// (`OQ`).
    OriginalQualityScores,
    /// (`OX`).
    OriginalUmiBarcodeSequence,
    /// (`PG`).
    Program,
    /// (`PQ`).
    TemplateLikelihood,
    /// (`PT`).
    PaddedReadAnnotations,
    /// (`PU`).
    PlatformUnit,
    /// (`Q2`).
    MateQualityScores,
    /// (`QT`).
    SampleBarcodeQualityScores,
    /// (`QX`).
    UmiQualityScores,
    /// (`R2`).
    MateSequence,
    /// (`RG`).
    ReadGroup,
    /// (`RT`).
    ReservedRt,
    /// (`RX`).
    UmiSequence,
    /// (`S2`).
    ReservedS2,
    /// (`SA`).
    OtherAlignments,
    /// (`SM`).
    TemplateMappingQuality,
    /// (`SQ`).
    ReservedSq,
    /// (`TC`).
    SegmentCount,
    /// (`TS`).
    TranscriptStrand,
    /// (`U2`).
    NextHitQualityScores,
    /// (`UQ`).
    SegmentLikelihood,
    /// Any other non-standard tag.
    Other(Other),
}

impl AsRef<[u8; LENGTH]> for Tag {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Self::MinMappingQuality => b"AM",
            Self::AlignmentScore => b"AS",
            Self::SampleBarcodeSequence => b"BC",
            Self::BaseAlignmentQualityOffsets => b"BQ",
            Self::OriginalUmiQualityScores => b"BZ",
            Self::CellBarcodeId => b"CB",
            Self::NextHitReferenceSequenceName => b"CC",
            Self::Cigar => b"CG",
            Self::ColorEditDistance => b"CM",
            Self::Comment => b"CO",
            Self::NextHitPosition => b"CP",
            Self::ColarQualityScores => b"CQ",
            Self::CellBarcodeSequence => b"CR",
            Self::ColorSequence => b"CS",
            Self::CompleteReadAnnotations => b"CT",
            Self::CellBarcodeQualityScores => b"CY",
            Self::NextHitSequence => b"E2",
            Self::SegmentIndex => b"FI",
            Self::SegmentSuffix => b"FS",
            Self::AlternativeSequence => b"FZ",
            Self::ReservedGc => b"GC",
            Self::ReservedGq => b"GQ",
            Self::ReservedGs => b"GS",
            Self::PerfectHitCount => b"HO",
            Self::OneDifferenceHitCount => b"H1",
            Self::TwoDifferenceHitCount => b"H2",
            Self::HitIndex => b"HI",
            Self::TotalHitCount => b"IH",
            Self::Library => b"LB",
            Self::MateCigar => b"MC",
            Self::MismatchedPositions => b"MD",
            Self::ReservedMf => b"MF",
            Self::UmiId => b"MI",
            Self::BaseModificationProbabilities => b"ML",
            Self::BaseModifications => b"MM",
            Self::MateMappingQuality => b"MQ",
            Self::AlignmentHitCount => b"NH",
            Self::EditDistance => b"NM",
            Self::OriginalAlignment => b"OA",
            Self::OriginalCigar => b"OC",
            Self::OriginalPosition => b"OP",
            Self::OriginalQualityScores => b"OQ",
            Self::OriginalUmiBarcodeSequence => b"OX",
            Self::Program => b"PG",
            Self::TemplateLikelihood => b"PQ",
            Self::PaddedReadAnnotations => b"PT",
            Self::PlatformUnit => b"PU",
            Self::MateQualityScores => b"Q2",
            Self::SampleBarcodeQualityScores => b"QT",
            Self::UmiQualityScores => b"QX",
            Self::MateSequence => b"R2",
            Self::ReadGroup => b"RG",
            Self::ReservedRt => b"RT",
            Self::UmiSequence => b"RX",
            Self::ReservedS2 => b"S2",
            Self::OtherAlignments => b"SA",
            Self::TemplateMappingQuality => b"SM",
            Self::ReservedSq => b"SQ",
            Self::SegmentCount => b"TC",
            Self::TranscriptStrand => b"TS",
            Self::NextHitQualityScores => b"U2",
            Self::SegmentLikelihood => b"UQ",
            Self::Other(Other(tag)) => tag,
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bytes = self.as_ref();
        f.write_char(char::from(bytes[0]))?;
        f.write_char(char::from(bytes[1]))?;
        Ok(())
    }
}

/// An error returned when a raw SAM record data field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The length is invalid.
    ///
    /// The tag length must be 2 characters.
    InvalidLength(usize),
    /// A character is invalid.
    InvalidCharacter(char),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLength(len) => {
                write!(f, "invalid length: expected {}, got {}", LENGTH, len)
            }
            Self::InvalidCharacter(c) => write!(f, "invalid character: {}", c),
        }
    }
}

// § 1.5 The alignment section: optional fields (2021-01-07)
impl TryFrom<[u8; LENGTH]> for Tag {
    type Error = ParseError;

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match &b {
            b"AM" => Ok(Self::MinMappingQuality),
            b"AS" => Ok(Self::AlignmentScore),
            b"BC" => Ok(Self::SampleBarcodeSequence),
            b"BQ" => Ok(Self::BaseAlignmentQualityOffsets),
            b"BZ" => Ok(Self::OriginalUmiQualityScores),
            b"CB" => Ok(Self::CellBarcodeId),
            b"CC" => Ok(Self::NextHitReferenceSequenceName),
            b"CG" => Ok(Self::Cigar),
            b"CM" => Ok(Self::ColorEditDistance),
            b"CO" => Ok(Self::Comment),
            b"CP" => Ok(Self::NextHitPosition),
            b"CQ" => Ok(Self::ColarQualityScores),
            b"CR" => Ok(Self::CellBarcodeSequence),
            b"CS" => Ok(Self::ColorSequence),
            b"CT" => Ok(Self::CompleteReadAnnotations),
            b"CY" => Ok(Self::CellBarcodeQualityScores),
            b"E2" => Ok(Self::NextHitSequence),
            b"FI" => Ok(Self::SegmentIndex),
            b"FS" => Ok(Self::SegmentSuffix),
            b"FZ" => Ok(Self::AlternativeSequence),
            b"GC" => Ok(Self::ReservedGc),
            b"GQ" => Ok(Self::ReservedGq),
            b"GS" => Ok(Self::ReservedGs),
            b"HO" => Ok(Self::PerfectHitCount),
            b"H1" => Ok(Self::OneDifferenceHitCount),
            b"H2" => Ok(Self::TwoDifferenceHitCount),
            b"HI" => Ok(Self::HitIndex),
            b"IH" => Ok(Self::TotalHitCount),
            b"LB" => Ok(Self::Library),
            b"MC" => Ok(Self::MateCigar),
            b"MD" => Ok(Self::MismatchedPositions),
            b"MF" => Ok(Self::ReservedMf),
            b"MI" => Ok(Self::UmiId),
            b"ML" => Ok(Self::BaseModificationProbabilities),
            b"MM" => Ok(Self::BaseModifications),
            b"MQ" => Ok(Self::MateMappingQuality),
            b"NH" => Ok(Self::AlignmentHitCount),
            b"NM" => Ok(Self::EditDistance),
            b"OA" => Ok(Self::OriginalAlignment),
            b"OC" => Ok(Self::OriginalCigar),
            b"OP" => Ok(Self::OriginalPosition),
            b"OQ" => Ok(Self::OriginalQualityScores),
            b"OX" => Ok(Self::OriginalUmiBarcodeSequence),
            b"PG" => Ok(Self::Program),
            b"PQ" => Ok(Self::TemplateLikelihood),
            b"PT" => Ok(Self::PaddedReadAnnotations),
            b"PU" => Ok(Self::PlatformUnit),
            b"Q2" => Ok(Self::MateQualityScores),
            b"QT" => Ok(Self::SampleBarcodeQualityScores),
            b"QX" => Ok(Self::UmiQualityScores),
            b"R2" => Ok(Self::MateSequence),
            b"RG" => Ok(Self::ReadGroup),
            b"RT" => Ok(Self::ReservedRt),
            b"RX" => Ok(Self::UmiSequence),
            b"S2" => Ok(Self::ReservedS2),
            b"SA" => Ok(Self::OtherAlignments),
            b"SM" => Ok(Self::TemplateMappingQuality),
            b"SQ" => Ok(Self::ReservedSq),
            b"TC" => Ok(Self::SegmentCount),
            b"TS" => Ok(Self::TranscriptStrand),
            b"U2" => Ok(Self::NextHitQualityScores),
            b"UQ" => Ok(Self::SegmentLikelihood),
            _ => {
                if !b[0].is_ascii_alphabetic() {
                    Err(ParseError::InvalidCharacter(char::from(b[0])))
                } else if !b[1].is_ascii_alphanumeric() {
                    Err(ParseError::InvalidCharacter(char::from(b[1])))
                } else {
                    Ok(Self::Other(Other(b)))
                }
            }
        }
    }
}

impl TryFrom<&[u8]> for Tag {
    type Error = ParseError;

    fn try_from(b: &[u8]) -> Result<Self, Self::Error> {
        if b.len() == LENGTH {
            Self::try_from([b[0], b[1]])
        } else {
            Err(ParseError::InvalidLength(b.len()))
        }
    }
}

impl FromStr for Tag {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Tag::try_from(s.as_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Tag::MinMappingQuality.to_string(), "AM");
        assert_eq!(Tag::AlignmentScore.to_string(), "AS");
        assert_eq!(Tag::SampleBarcodeSequence.to_string(), "BC");
        assert_eq!(Tag::BaseAlignmentQualityOffsets.to_string(), "BQ");
        assert_eq!(Tag::OriginalUmiQualityScores.to_string(), "BZ");
        assert_eq!(Tag::CellBarcodeId.to_string(), "CB");
        assert_eq!(Tag::NextHitReferenceSequenceName.to_string(), "CC");
        assert_eq!(Tag::Cigar.to_string(), "CG");
        assert_eq!(Tag::ColorEditDistance.to_string(), "CM");
        assert_eq!(Tag::Comment.to_string(), "CO");
        assert_eq!(Tag::NextHitPosition.to_string(), "CP");
        assert_eq!(Tag::ColarQualityScores.to_string(), "CQ");
        assert_eq!(Tag::CellBarcodeSequence.to_string(), "CR");
        assert_eq!(Tag::ColorSequence.to_string(), "CS");
        assert_eq!(Tag::CompleteReadAnnotations.to_string(), "CT");
        assert_eq!(Tag::CellBarcodeQualityScores.to_string(), "CY");
        assert_eq!(Tag::NextHitSequence.to_string(), "E2");
        assert_eq!(Tag::SegmentIndex.to_string(), "FI");
        assert_eq!(Tag::SegmentSuffix.to_string(), "FS");
        assert_eq!(Tag::AlternativeSequence.to_string(), "FZ");
        assert_eq!(Tag::ReservedGc.to_string(), "GC");
        assert_eq!(Tag::ReservedGq.to_string(), "GQ");
        assert_eq!(Tag::ReservedGs.to_string(), "GS");
        assert_eq!(Tag::PerfectHitCount.to_string(), "HO");
        assert_eq!(Tag::OneDifferenceHitCount.to_string(), "H1");
        assert_eq!(Tag::TwoDifferenceHitCount.to_string(), "H2");
        assert_eq!(Tag::HitIndex.to_string(), "HI");
        assert_eq!(Tag::TotalHitCount.to_string(), "IH");
        assert_eq!(Tag::Library.to_string(), "LB");
        assert_eq!(Tag::MateCigar.to_string(), "MC");
        assert_eq!(Tag::MismatchedPositions.to_string(), "MD");
        assert_eq!(Tag::ReservedMf.to_string(), "MF");
        assert_eq!(Tag::UmiId.to_string(), "MI");
        assert_eq!(Tag::BaseModificationProbabilities.to_string(), "ML");
        assert_eq!(Tag::BaseModifications.to_string(), "MM");
        assert_eq!(Tag::MateMappingQuality.to_string(), "MQ");
        assert_eq!(Tag::AlignmentHitCount.to_string(), "NH");
        assert_eq!(Tag::EditDistance.to_string(), "NM");
        assert_eq!(Tag::OriginalAlignment.to_string(), "OA");
        assert_eq!(Tag::OriginalCigar.to_string(), "OC");
        assert_eq!(Tag::OriginalPosition.to_string(), "OP");
        assert_eq!(Tag::OriginalQualityScores.to_string(), "OQ");
        assert_eq!(Tag::OriginalUmiBarcodeSequence.to_string(), "OX");
        assert_eq!(Tag::Program.to_string(), "PG");
        assert_eq!(Tag::TemplateLikelihood.to_string(), "PQ");
        assert_eq!(Tag::PaddedReadAnnotations.to_string(), "PT");
        assert_eq!(Tag::PlatformUnit.to_string(), "PU");
        assert_eq!(Tag::MateQualityScores.to_string(), "Q2");
        assert_eq!(Tag::SampleBarcodeQualityScores.to_string(), "QT");
        assert_eq!(Tag::UmiQualityScores.to_string(), "QX");
        assert_eq!(Tag::MateSequence.to_string(), "R2");
        assert_eq!(Tag::ReadGroup.to_string(), "RG");
        assert_eq!(Tag::ReservedRt.to_string(), "RT");
        assert_eq!(Tag::UmiSequence.to_string(), "RX");
        assert_eq!(Tag::ReservedS2.to_string(), "S2");
        assert_eq!(Tag::OtherAlignments.to_string(), "SA");
        assert_eq!(Tag::TemplateMappingQuality.to_string(), "SM");
        assert_eq!(Tag::ReservedSq.to_string(), "SQ");
        assert_eq!(Tag::SegmentCount.to_string(), "TC");
        assert_eq!(Tag::TranscriptStrand.to_string(), "TS");
        assert_eq!(Tag::NextHitQualityScores.to_string(), "U2");
        assert_eq!(Tag::SegmentLikelihood.to_string(), "UQ");
        assert_eq!(Tag::Other(Other([b'Z', b'N'])).to_string(), "ZN");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("AM".parse(), Ok(Tag::MinMappingQuality));
        assert_eq!("AS".parse(), Ok(Tag::AlignmentScore));
        assert_eq!("BC".parse(), Ok(Tag::SampleBarcodeSequence));
        assert_eq!("BQ".parse(), Ok(Tag::BaseAlignmentQualityOffsets));
        assert_eq!("BZ".parse(), Ok(Tag::OriginalUmiQualityScores));
        assert_eq!("CB".parse(), Ok(Tag::CellBarcodeId));
        assert_eq!("CC".parse(), Ok(Tag::NextHitReferenceSequenceName));
        assert_eq!("CG".parse(), Ok(Tag::Cigar));
        assert_eq!("CM".parse(), Ok(Tag::ColorEditDistance));
        assert_eq!("CO".parse(), Ok(Tag::Comment));
        assert_eq!("CP".parse(), Ok(Tag::NextHitPosition));
        assert_eq!("CQ".parse(), Ok(Tag::ColarQualityScores));
        assert_eq!("CR".parse(), Ok(Tag::CellBarcodeSequence));
        assert_eq!("CS".parse(), Ok(Tag::ColorSequence));
        assert_eq!("CT".parse(), Ok(Tag::CompleteReadAnnotations));
        assert_eq!("CY".parse(), Ok(Tag::CellBarcodeQualityScores));
        assert_eq!("E2".parse(), Ok(Tag::NextHitSequence));
        assert_eq!("FI".parse(), Ok(Tag::SegmentIndex));
        assert_eq!("FS".parse(), Ok(Tag::SegmentSuffix));
        assert_eq!("FZ".parse(), Ok(Tag::AlternativeSequence));
        assert_eq!("GC".parse(), Ok(Tag::ReservedGc));
        assert_eq!("GQ".parse(), Ok(Tag::ReservedGq));
        assert_eq!("GS".parse(), Ok(Tag::ReservedGs));
        assert_eq!("HO".parse(), Ok(Tag::PerfectHitCount));
        assert_eq!("H1".parse(), Ok(Tag::OneDifferenceHitCount));
        assert_eq!("H2".parse(), Ok(Tag::TwoDifferenceHitCount));
        assert_eq!("HI".parse(), Ok(Tag::HitIndex));
        assert_eq!("IH".parse(), Ok(Tag::TotalHitCount));
        assert_eq!("LB".parse(), Ok(Tag::Library));
        assert_eq!("MC".parse(), Ok(Tag::MateCigar));
        assert_eq!("MD".parse(), Ok(Tag::MismatchedPositions));
        assert_eq!("MF".parse(), Ok(Tag::ReservedMf));
        assert_eq!("MI".parse(), Ok(Tag::UmiId));
        assert_eq!("ML".parse(), Ok(Tag::BaseModificationProbabilities));
        assert_eq!("MM".parse(), Ok(Tag::BaseModifications));
        assert_eq!("MQ".parse(), Ok(Tag::MateMappingQuality));
        assert_eq!("NH".parse(), Ok(Tag::AlignmentHitCount));
        assert_eq!("NM".parse(), Ok(Tag::EditDistance));
        assert_eq!("OA".parse(), Ok(Tag::OriginalAlignment));
        assert_eq!("OC".parse(), Ok(Tag::OriginalCigar));
        assert_eq!("OP".parse(), Ok(Tag::OriginalPosition));
        assert_eq!("OQ".parse(), Ok(Tag::OriginalQualityScores));
        assert_eq!("OX".parse(), Ok(Tag::OriginalUmiBarcodeSequence));
        assert_eq!("PG".parse(), Ok(Tag::Program));
        assert_eq!("PQ".parse(), Ok(Tag::TemplateLikelihood));
        assert_eq!("PT".parse(), Ok(Tag::PaddedReadAnnotations));
        assert_eq!("PU".parse(), Ok(Tag::PlatformUnit));
        assert_eq!("Q2".parse(), Ok(Tag::MateQualityScores));
        assert_eq!("QT".parse(), Ok(Tag::SampleBarcodeQualityScores));
        assert_eq!("QX".parse(), Ok(Tag::UmiQualityScores));
        assert_eq!("R2".parse(), Ok(Tag::MateSequence));
        assert_eq!("RG".parse(), Ok(Tag::ReadGroup));
        assert_eq!("RT".parse(), Ok(Tag::ReservedRt));
        assert_eq!("RX".parse(), Ok(Tag::UmiSequence));
        assert_eq!("S2".parse(), Ok(Tag::ReservedS2));
        assert_eq!("SA".parse(), Ok(Tag::OtherAlignments));
        assert_eq!("SM".parse(), Ok(Tag::TemplateMappingQuality));
        assert_eq!("SQ".parse(), Ok(Tag::ReservedSq));
        assert_eq!("TC".parse(), Ok(Tag::SegmentCount));
        assert_eq!("TS".parse(), Ok(Tag::TranscriptStrand));
        assert_eq!("U2".parse(), Ok(Tag::NextHitQualityScores));
        assert_eq!("UQ".parse(), Ok(Tag::SegmentLikelihood));
        assert_eq!("ZN".parse(), Ok(Tag::Other(Other([b'Z', b'N']))));

        assert_eq!("".parse::<Tag>(), Err(ParseError::InvalidLength(0)));
        assert_eq!("R".parse::<Tag>(), Err(ParseError::InvalidLength(1)));
        assert_eq!("RGP".parse::<Tag>(), Err(ParseError::InvalidLength(3)));
        assert_eq!("1G".parse::<Tag>(), Err(ParseError::InvalidCharacter('1')));
        assert_eq!("R_".parse::<Tag>(), Err(ParseError::InvalidCharacter('_')));
    }
}
