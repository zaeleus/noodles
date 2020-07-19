use std::{convert::TryFrom, error, fmt};

use crate::num::Itf8;

#[derive(Debug, Eq, PartialEq)]
pub struct TryFromItf8Error(Itf8);

impl fmt::Display for TryFromItf8Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid encoding kind: expected 0..=9, got {}", self.0)
    }
}

impl error::Error for TryFromItf8Error {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    Null,
    External,
    Golomb,
    Huffman,
    ByteArrayLen,
    ByteArrayStop,
    Beta,
    Subexp,
    GolombRice,
    Gamma,
}

impl TryFrom<Itf8> for Kind {
    type Error = TryFromItf8Error;

    fn try_from(id: Itf8) -> Result<Self, Self::Error> {
        match id {
            0 => Ok(Self::Null),
            1 => Ok(Self::External),
            2 => Ok(Self::Golomb),
            3 => Ok(Self::Huffman),
            4 => Ok(Self::ByteArrayLen),
            5 => Ok(Self::ByteArrayStop),
            6 => Ok(Self::Beta),
            7 => Ok(Self::Subexp),
            8 => Ok(Self::GolombRice),
            9 => Ok(Self::Gamma),
            _ => Err(TryFromItf8Error(id)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from() {
        assert_eq!(Kind::try_from(0), Ok(Kind::Null));
        assert_eq!(Kind::try_from(1), Ok(Kind::External));
        assert_eq!(Kind::try_from(2), Ok(Kind::Golomb));
        assert_eq!(Kind::try_from(3), Ok(Kind::Huffman));
        assert_eq!(Kind::try_from(4), Ok(Kind::ByteArrayLen));
        assert_eq!(Kind::try_from(5), Ok(Kind::ByteArrayStop));
        assert_eq!(Kind::try_from(6), Ok(Kind::Beta));
        assert_eq!(Kind::try_from(7), Ok(Kind::Subexp));
        assert_eq!(Kind::try_from(8), Ok(Kind::GolombRice));
        assert_eq!(Kind::try_from(9), Ok(Kind::Gamma));
        assert_eq!(Kind::try_from(10), Err(TryFromItf8Error(10)));
    }
}
