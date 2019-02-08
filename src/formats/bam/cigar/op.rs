use std::fmt;

use byteorder::{ByteOrder, LittleEndian};

#[derive(Debug, Eq, PartialEq)]
pub enum Op {
    Match(u32),
    Insertion(u32),
    Deletion(u32),
    Skip(u32),
    SoftClip(u32),
    HardClip(u32),
    Pad(u32),
    SeqMatch(u32),
    SeqMismatch(u32),
}

impl Op {
    pub fn from_bytes(bytes: &[u8]) -> Op {
        Op::from_u32(LittleEndian::read_u32(bytes))
    }

    pub fn from_u32(u: u32) -> Op {
        let len = u >> 4;

        match u & 0x0f {
            0 => Op::Match(len),
            1 => Op::Insertion(len),
            2 => Op::Deletion(len),
            3 => Op::Skip(len),
            4 => Op::SoftClip(len),
            5 => Op::HardClip(len),
            6 => Op::Pad(len),
            7 => Op::SeqMatch(len),
            8 => Op::SeqMismatch(len),
            _ => panic!("invalid CIGAR op"),
        }
    }

    pub fn len(&self) -> u32 {
        match *self {
            Op::Match(len) => len,
            Op::Insertion(len) => len,
            Op::Deletion(len) => len,
            Op::Skip(len) => len,
            Op::SoftClip(len) => len,
            Op::HardClip(len) => len,
            Op::Pad(len) => len,
            Op::SeqMatch(len) => len,
            Op::SeqMismatch(len) => len,
        }
    }

    pub fn op(&self) -> char {
        match *self {
            Op::Match(_) => 'M',
            Op::Insertion(_) => 'I',
            Op::Deletion(_) => 'D',
            Op::Skip(_) => 'N',
            Op::SoftClip(_) => 'S',
            Op::HardClip(_) => 'H',
            Op::Pad(_) => 'P',
            Op::SeqMatch(_) => '=',
            Op::SeqMismatch(_) => 'X',
        }
    }
}

impl fmt::Display for Op {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{}", self.len(), self.op())
    }
}

#[cfg(test)]
mod tests {
    use super::Op;

    #[test]
    fn test_from_u32() {
        assert_eq!(Op::from_u32(1 << 4 | 0), Op::Match(1));
        assert_eq!(Op::from_u32(2 << 4 | 1), Op::Insertion(2));
        assert_eq!(Op::from_u32(3 << 4 | 2), Op::Deletion(3));
        assert_eq!(Op::from_u32(4 << 4 | 3), Op::Skip(4));
        assert_eq!(Op::from_u32(5 << 4 | 4), Op::SoftClip(5));
        assert_eq!(Op::from_u32(6 << 4 | 5), Op::HardClip(6));
        assert_eq!(Op::from_u32(7 << 4 | 6), Op::Pad(7));
        assert_eq!(Op::from_u32(8 << 4 | 7), Op::SeqMatch(8));
        assert_eq!(Op::from_u32(9 << 4 | 8), Op::SeqMismatch(9));
    }

    #[test]
    #[should_panic]
    fn test_from_u32_with_invalid_op() {
        Op::from_u32(1 << 4 | 9);
    }

    #[test]
    fn test_len() {
        assert_eq!(Op::Match(1).len(), 1);
        assert_eq!(Op::Insertion(2).len(), 2);
        assert_eq!(Op::Deletion(3).len(), 3);
        assert_eq!(Op::Skip(4).len(), 4);
        assert_eq!(Op::SoftClip(5).len(), 5);
        assert_eq!(Op::HardClip(6).len(), 6);
        assert_eq!(Op::Pad(7).len(), 7);
        assert_eq!(Op::SeqMatch(8).len(), 8);
        assert_eq!(Op::SeqMismatch(9).len(), 9);
    }

    #[test]
    fn test_op() {
        assert_eq!(Op::Match(1).op(), 'M');
        assert_eq!(Op::Insertion(1).op(), 'I');
        assert_eq!(Op::Deletion(1).op(), 'D');
        assert_eq!(Op::Skip(1).op(), 'N');
        assert_eq!(Op::SoftClip(1).op(), 'S');
        assert_eq!(Op::HardClip(1).op(), 'H');
        assert_eq!(Op::Pad(1).op(), 'P');
        assert_eq!(Op::SeqMatch(1).op(), '=');
        assert_eq!(Op::SeqMismatch(1).op(), 'X');
    }

    #[test]
    fn test_fmt() {
        assert_eq!(format!("{}", Op::Match(32)), "32M");
        assert_eq!(format!("{}", Op::Deletion(7)), "7D");
        assert_eq!(format!("{}", Op::Skip(11)), "11N");
        assert_eq!(format!("{}", Op::Pad(188)), "188P");
    }
}
