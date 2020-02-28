pub use self::kind::Kind;

pub mod kind;

use std::{convert::TryFrom, error, fmt};

use byteorder::{ByteOrder, LittleEndian};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Op {
    kind: Kind,
    len: u32,
}

impl Op {
    pub fn new(kind: Kind, len: u32) -> Self {
        Self { kind, len }
    }

    pub fn kind(self) -> Kind {
        self.kind
    }

    pub fn len(self) -> u32 {
        self.len
    }

    pub fn is_empty(self) -> bool {
        self.len == 0
    }
}

#[derive(Debug, Eq, PartialEq)]
pub struct TryFromUintError(u32);

impl fmt::Display for TryFromUintError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid op: expected 0..=8, got {}", self.0)
    }
}

impl error::Error for TryFromUintError {}

impl TryFrom<u32> for Op {
    type Error = TryFromUintError;

    fn try_from(u: u32) -> Result<Self, Self::Error> {
        let len = u >> 4;

        let kind = match u & 0x0f {
            0 => Kind::Match,
            1 => Kind::Insertion,
            2 => Kind::Deletion,
            3 => Kind::Skip,
            4 => Kind::SoftClip,
            5 => Kind::HardClip,
            6 => Kind::Pad,
            7 => Kind::SeqMatch,
            8 => Kind::SeqMismatch,
            n => return Err(TryFromUintError(n)),
        };

        Ok(Self::new(kind, len))
    }
}

impl TryFrom<&[u8]> for Op {
    type Error = ();

    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        let u = LittleEndian::read_u32(bytes);
        Self::try_from(u).map_err(|_| ())
    }
}

impl fmt::Display for Op {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}", self.len(), self.kind())
    }
}

impl From<Op> for u32 {
    fn from(op: Op) -> u32 {
        let i = op.kind() as u32;
        op.len() << 4 | i
    }
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use crate::cigar::Kind;

    use super::*;

    #[test]
    fn test_new() {
        let op = Op::new(Kind::Match, 13);
        assert_eq!(op.kind, Kind::Match);
        assert_eq!(op.len, 13);
    }

    #[test]
    fn test_kind() {
        let op = Op::new(Kind::SoftClip, 13);
        assert_eq!(op.kind(), Kind::SoftClip);
    }

    #[test]
    fn test_len() {
        let op = Op::new(Kind::Match, 1);
        assert_eq!(op.len(), 1);
    }

    #[test]
    fn test_is_empty() {
        let op = Op::new(Kind::Match, 0);
        assert!(op.is_empty());

        let op = Op::new(Kind::Match, 1);
        assert!(!op.is_empty());
    }

    #[test]
    fn test_try_from_u32() {
        assert_eq!(Op::try_from(1 << 4 | 0), Ok(Op::new(Kind::Match, 1)));
        assert_eq!(Op::try_from(2 << 4 | 1), Ok(Op::new(Kind::Insertion, 2)));
        assert_eq!(Op::try_from(3 << 4 | 2), Ok(Op::new(Kind::Deletion, 3)));
        assert_eq!(Op::try_from(4 << 4 | 3), Ok(Op::new(Kind::Skip, 4)));
        assert_eq!(Op::try_from(5 << 4 | 4), Ok(Op::new(Kind::SoftClip, 5)));
        assert_eq!(Op::try_from(6 << 4 | 5), Ok(Op::new(Kind::HardClip, 6)));
        assert_eq!(Op::try_from(7 << 4 | 6), Ok(Op::new(Kind::Pad, 7)));
        assert_eq!(Op::try_from(8 << 4 | 7), Ok(Op::new(Kind::SeqMatch, 8)));
        assert_eq!(Op::try_from(9 << 4 | 8), Ok(Op::new(Kind::SeqMismatch, 9)));
        assert_eq!(Op::try_from(10 << 4 | 9), Err(TryFromUintError(9)));
    }

    #[test]
    fn test_try_from_u8_slice() {
        let bytes = [0x40, 0x02, 0x00, 0x00];
        assert_eq!(Op::try_from(&bytes[..]), Ok(Op::new(Kind::Match, 36)));
    }

    #[test]
    fn test_fmt() {
        assert_eq!(format!("{}", Op::new(Kind::Match, 32)), "32M");
        assert_eq!(format!("{}", Op::new(Kind::Deletion, 7)), "7D");
        assert_eq!(format!("{}", Op::new(Kind::Skip, 11)), "11N");
        assert_eq!(format!("{}", Op::new(Kind::Pad, 188)), "188P");
    }

    #[test]
    fn test_from_op_for_u32() {
        assert_eq!(u32::from(Op::new(Kind::Match, 1)), 1 << 4 | 0);
        assert_eq!(u32::from(Op::new(Kind::Insertion, 2)), 2 << 4 | 1);
        assert_eq!(u32::from(Op::new(Kind::Deletion, 3)), 3 << 4 | 2);
        assert_eq!(u32::from(Op::new(Kind::Skip, 4)), 4 << 4 | 3);
        assert_eq!(u32::from(Op::new(Kind::SoftClip, 5)), 5 << 4 | 4);
        assert_eq!(u32::from(Op::new(Kind::HardClip, 6)), 6 << 4 | 5);
        assert_eq!(u32::from(Op::new(Kind::Pad, 7)), 7 << 4 | 6);
        assert_eq!(u32::from(Op::new(Kind::SeqMatch, 8)), 8 << 4 | 7);
        assert_eq!(u32::from(Op::new(Kind::SeqMismatch, 9)), 9 << 4 | 8);
    }
}
