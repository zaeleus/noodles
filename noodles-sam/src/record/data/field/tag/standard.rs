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

impl Standard {
    pub(super) const fn new(buf: [u8; LENGTH]) -> Option<Self> {
        match &buf {
            b"AM" => Some(Self::MinMappingQuality),
            b"AS" => Some(Self::AlignmentScore),
            b"BC" => Some(Self::SampleBarcodeSequence),
            b"BQ" => Some(Self::BaseAlignmentQualityOffsets),
            b"BZ" => Some(Self::OriginalUmiQualityScores),
            b"CB" => Some(Self::CellBarcodeId),
            b"CC" => Some(Self::NextHitReferenceSequenceName),
            b"CG" => Some(Self::Cigar),
            b"CM" => Some(Self::ColorEditDistance),
            b"CO" => Some(Self::Comment),
            b"CP" => Some(Self::NextHitPosition),
            b"CQ" => Some(Self::ColarQualityScores),
            b"CR" => Some(Self::CellBarcodeSequence),
            b"CS" => Some(Self::ColorSequence),
            b"CT" => Some(Self::CompleteReadAnnotations),
            b"CY" => Some(Self::CellBarcodeQualityScores),
            b"E2" => Some(Self::NextHitSequence),
            b"FI" => Some(Self::SegmentIndex),
            b"FS" => Some(Self::SegmentSuffix),
            b"FZ" => Some(Self::AlternativeSequence),
            b"GC" => Some(Self::ReservedGc),
            b"GQ" => Some(Self::ReservedGq),
            b"GS" => Some(Self::ReservedGs),
            b"HO" => Some(Self::PerfectHitCount),
            b"H1" => Some(Self::OneDifferenceHitCount),
            b"H2" => Some(Self::TwoDifferenceHitCount),
            b"HI" => Some(Self::HitIndex),
            b"IH" => Some(Self::TotalHitCount),
            b"LB" => Some(Self::Library),
            b"MC" => Some(Self::MateCigar),
            b"MD" => Some(Self::MismatchedPositions),
            b"MF" => Some(Self::ReservedMf),
            b"MI" => Some(Self::UmiId),
            b"ML" => Some(Self::BaseModificationProbabilities),
            b"MM" => Some(Self::BaseModifications),
            b"MN" => Some(Self::BaseModificationSequenceLength),
            b"MQ" => Some(Self::MateMappingQuality),
            b"NH" => Some(Self::AlignmentHitCount),
            b"NM" => Some(Self::EditDistance),
            b"OA" => Some(Self::OriginalAlignment),
            b"OC" => Some(Self::OriginalCigar),
            b"OP" => Some(Self::OriginalPosition),
            b"OQ" => Some(Self::OriginalQualityScores),
            b"OX" => Some(Self::OriginalUmiBarcodeSequence),
            b"PG" => Some(Self::Program),
            b"PQ" => Some(Self::TemplateLikelihood),
            b"PT" => Some(Self::PaddedReadAnnotations),
            b"PU" => Some(Self::PlatformUnit),
            b"Q2" => Some(Self::MateQualityScores),
            b"QT" => Some(Self::SampleBarcodeQualityScores),
            b"QX" => Some(Self::UmiQualityScores),
            b"R2" => Some(Self::MateSequence),
            b"RG" => Some(Self::ReadGroup),
            b"RT" => Some(Self::ReservedRt),
            b"RX" => Some(Self::UmiSequence),
            b"S2" => Some(Self::ReservedS2),
            b"SA" => Some(Self::OtherAlignments),
            b"SM" => Some(Self::TemplateMappingQuality),
            b"SQ" => Some(Self::ReservedSq),
            b"TC" => Some(Self::SegmentCount),
            b"TS" => Some(Self::TranscriptStrand),
            b"U2" => Some(Self::NextHitQualityScores),
            b"UQ" => Some(Self::SegmentLikelihood),
            _ => None,
        }
    }
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
