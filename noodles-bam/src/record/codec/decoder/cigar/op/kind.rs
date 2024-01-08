use std::{error, fmt};

use noodles_sam::alignment::record::cigar::op::Kind;

/// An error when a raw BAM record CIGAR op kind fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// The input is invalid.
    Invalid { actual: u8 },
}

impl error::Error for DecodeError {}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid { actual } => write!(f, "invalid input: expected 0..=8, got {actual}"),
        }
    }
}

pub(super) fn decode_kind(n: u32) -> Result<Kind, DecodeError> {
    match n & 0x0f {
        0 => Ok(Kind::Match),
        1 => Ok(Kind::Insertion),
        2 => Ok(Kind::Deletion),
        3 => Ok(Kind::Skip),
        4 => Ok(Kind::SoftClip),
        5 => Ok(Kind::HardClip),
        6 => Ok(Kind::Pad),
        7 => Ok(Kind::SequenceMatch),
        8 => Ok(Kind::SequenceMismatch),
        n => Err(DecodeError::Invalid {
            // SAFETY: `n` is <= 15.
            actual: n as u8,
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_kind() {
        assert_eq!(decode_kind(0x00), Ok(Kind::Match));
        assert_eq!(decode_kind(0x01), Ok(Kind::Insertion));
        assert_eq!(decode_kind(0x02), Ok(Kind::Deletion));
        assert_eq!(decode_kind(0x03), Ok(Kind::Skip));
        assert_eq!(decode_kind(0x04), Ok(Kind::SoftClip));
        assert_eq!(decode_kind(0x05), Ok(Kind::HardClip));
        assert_eq!(decode_kind(0x06), Ok(Kind::Pad));
        assert_eq!(decode_kind(0x07), Ok(Kind::SequenceMatch));
        assert_eq!(decode_kind(0x08), Ok(Kind::SequenceMismatch));

        assert_eq!(decode_kind(0x09), Err(DecodeError::Invalid { actual: 9 }));
    }
}
