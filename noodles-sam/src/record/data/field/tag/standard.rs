use super::LENGTH;

/// Standard tags are defined in "Sequence Alignment/Map Optional Fields Specification"
/// (2020-05-29).
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum Standard {
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
    /// (`MN`).
    BaseModificationSequenceLength,
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
}

impl AsRef<[u8; LENGTH]> for Standard {
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
            Self::BaseModificationSequenceLength => b"MN",
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
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Standard {
    type Error = ();

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
            b"MN" => Ok(Self::BaseModificationSequenceLength),
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
            _ => Err(()),
        }
    }
}
