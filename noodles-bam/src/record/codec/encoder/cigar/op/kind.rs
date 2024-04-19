use noodles_sam::alignment::record::cigar::op::Kind;

pub(super) fn encode_kind(kind: Kind) -> u32 {
    match kind {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_kind() {
        assert_eq!(encode_kind(Kind::Match), 0);
        assert_eq!(encode_kind(Kind::Insertion), 1);
        assert_eq!(encode_kind(Kind::Deletion), 2);
        assert_eq!(encode_kind(Kind::Skip), 3);
        assert_eq!(encode_kind(Kind::SoftClip), 4);
        assert_eq!(encode_kind(Kind::HardClip), 5);
        assert_eq!(encode_kind(Kind::Pad), 6);
        assert_eq!(encode_kind(Kind::SequenceMatch), 7);
        assert_eq!(encode_kind(Kind::SequenceMismatch), 8);
    }
}
