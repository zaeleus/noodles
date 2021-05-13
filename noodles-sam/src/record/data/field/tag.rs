//! SAM record data field tag.

use std::{error, fmt, str::FromStr};

const LEN: usize = 2;

/// A SAM record data field tag.
///
/// Standard tags are defined in "Sequence Alignment/Map Optional Fields Specification"
/// (2020-05-29).
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
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
    Other(String),
}

impl AsRef<str> for Tag {
    fn as_ref(&self) -> &str {
        match self {
            Self::MinMappingQuality => "AM",
            Self::AlignmentScore => "AS",
            Self::SampleBarcodeSequence => "BC",
            Self::BaseAlignmentQualityOffsets => "BQ",
            Self::OriginalUmiQualityScores => "BZ",
            Self::CellBarcodeId => "CB",
            Self::NextHitReferenceSequenceName => "CC",
            Self::Cigar => "CG",
            Self::ColorEditDistance => "CM",
            Self::Comment => "CO",
            Self::NextHitPosition => "CP",
            Self::ColarQualityScores => "CQ",
            Self::CellBarcodeSequence => "CR",
            Self::ColorSequence => "CS",
            Self::CompleteReadAnnotations => "CT",
            Self::CellBarcodeQualityScores => "CY",
            Self::NextHitSequence => "E2",
            Self::SegmentIndex => "FI",
            Self::SegmentSuffix => "FS",
            Self::AlternativeSequence => "FZ",
            Self::ReservedGc => "GC",
            Self::ReservedGq => "GQ",
            Self::ReservedGs => "GS",
            Self::PerfectHitCount => "HO",
            Self::OneDifferenceHitCount => "H1",
            Self::TwoDifferenceHitCount => "H2",
            Self::HitIndex => "HI",
            Self::TotalHitCount => "IH",
            Self::Library => "LB",
            Self::MateCigar => "MC",
            Self::MismatchedPositions => "MD",
            Self::ReservedMf => "MF",
            Self::UmiId => "MI",
            Self::MateMappingQuality => "MQ",
            Self::AlignmentHitCount => "NH",
            Self::EditDistance => "NM",
            Self::OriginalAlignment => "OA",
            Self::OriginalCigar => "OC",
            Self::OriginalPosition => "OP",
            Self::OriginalQualityScores => "OQ",
            Self::OriginalUmiBarcodeSequence => "OX",
            Self::Program => "PG",
            Self::TemplateLikelihood => "PQ",
            Self::PaddedReadAnnotations => "PT",
            Self::PlatformUnit => "PU",
            Self::MateQualityScores => "Q2",
            Self::SampleBarcodeQualityScores => "QT",
            Self::UmiQualityScores => "QX",
            Self::MateSequence => "R2",
            Self::ReadGroup => "RG",
            Self::ReservedRt => "RT",
            Self::UmiSequence => "RX",
            Self::ReservedS2 => "S2",
            Self::OtherAlignments => "SA",
            Self::TemplateMappingQuality => "SM",
            Self::ReservedSq => "SQ",
            Self::SegmentCount => "TC",
            Self::TranscriptStrand => "TS",
            Self::NextHitQualityScores => "U2",
            Self::SegmentLikelihood => "UQ",
            Self::Other(tag) => tag,
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw SAM record data field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl FromStr for Tag {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "AM" => Ok(Self::MinMappingQuality),
            "AS" => Ok(Self::AlignmentScore),
            "BC" => Ok(Self::SampleBarcodeSequence),
            "BQ" => Ok(Self::BaseAlignmentQualityOffsets),
            "BZ" => Ok(Self::OriginalUmiQualityScores),
            "CB" => Ok(Self::CellBarcodeId),
            "CC" => Ok(Self::NextHitReferenceSequenceName),
            "CG" => Ok(Self::Cigar),
            "CM" => Ok(Self::ColorEditDistance),
            "CO" => Ok(Self::Comment),
            "CP" => Ok(Self::NextHitPosition),
            "CQ" => Ok(Self::ColarQualityScores),
            "CR" => Ok(Self::CellBarcodeSequence),
            "CS" => Ok(Self::ColorSequence),
            "CT" => Ok(Self::CompleteReadAnnotations),
            "CY" => Ok(Self::CellBarcodeQualityScores),
            "E2" => Ok(Self::NextHitSequence),
            "FI" => Ok(Self::SegmentIndex),
            "FS" => Ok(Self::SegmentSuffix),
            "FZ" => Ok(Self::AlternativeSequence),
            "GC" => Ok(Self::ReservedGc),
            "GQ" => Ok(Self::ReservedGq),
            "GS" => Ok(Self::ReservedGs),
            "HO" => Ok(Self::PerfectHitCount),
            "H1" => Ok(Self::OneDifferenceHitCount),
            "H2" => Ok(Self::TwoDifferenceHitCount),
            "HI" => Ok(Self::HitIndex),
            "IH" => Ok(Self::TotalHitCount),
            "LB" => Ok(Self::Library),
            "MC" => Ok(Self::MateCigar),
            "MD" => Ok(Self::MismatchedPositions),
            "MF" => Ok(Self::ReservedMf),
            "MI" => Ok(Self::UmiId),
            "MQ" => Ok(Self::MateMappingQuality),
            "NH" => Ok(Self::AlignmentHitCount),
            "NM" => Ok(Self::EditDistance),
            "OA" => Ok(Self::OriginalAlignment),
            "OC" => Ok(Self::OriginalCigar),
            "OP" => Ok(Self::OriginalPosition),
            "OQ" => Ok(Self::OriginalQualityScores),
            "OX" => Ok(Self::OriginalUmiBarcodeSequence),
            "PG" => Ok(Self::Program),
            "PQ" => Ok(Self::TemplateLikelihood),
            "PT" => Ok(Self::PaddedReadAnnotations),
            "PU" => Ok(Self::PlatformUnit),
            "Q2" => Ok(Self::MateQualityScores),
            "QT" => Ok(Self::SampleBarcodeQualityScores),
            "QX" => Ok(Self::UmiQualityScores),
            "R2" => Ok(Self::MateSequence),
            "RG" => Ok(Self::ReadGroup),
            "RT" => Ok(Self::ReservedRt),
            "RX" => Ok(Self::UmiSequence),
            "S2" => Ok(Self::ReservedS2),
            "SA" => Ok(Self::OtherAlignments),
            "SM" => Ok(Self::TemplateMappingQuality),
            "SQ" => Ok(Self::ReservedSq),
            "TC" => Ok(Self::SegmentCount),
            "TS" => Ok(Self::TranscriptStrand),
            "U2" => Ok(Self::NextHitQualityScores),
            "UQ" => Ok(Self::SegmentLikelihood),
            _ => {
                if s.len() == LEN {
                    Ok(Self::Other(s.into()))
                } else {
                    Err(ParseError(s.into()))
                }
            }
        }
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
        assert_eq!(Tag::Other(String::from("ZN")).to_string(), "ZN");
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
        assert_eq!("ZN".parse(), Ok(Tag::Other(String::from("ZN"))));

        assert_eq!("".parse::<Tag>(), Err(ParseError(String::from(""))));
        assert_eq!("R".parse::<Tag>(), Err(ParseError(String::from("R"))));
        assert_eq!("RGP".parse::<Tag>(), Err(ParseError(String::from("RGP"))));
        assert_eq!("RGRP".parse::<Tag>(), Err(ParseError(String::from("RGRP"))));
    }
}
