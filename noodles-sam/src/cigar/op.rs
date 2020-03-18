mod kind;

use std::{error, fmt, num, str::FromStr};

pub use self::kind::Kind;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Op {
    kind: Kind,
    len: u32,
}

#[allow(clippy::len_without_is_empty)]
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
}

#[derive(Debug)]
pub enum ParseError {
    InvalidLength(num::ParseIntError),
    InvalidKind(kind::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLength(e) => write!(f, "invalid op length: {}", e),
            Self::InvalidKind(e) => write!(f, "invalid op kind: {}", e),
        }
    }
}

impl FromStr for Op {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (raw_len, raw_kind) = s.split_at(s.len() - 1);
        let len = raw_len.parse().map_err(ParseError::InvalidLength)?;
        let kind = raw_kind.parse().map_err(ParseError::InvalidKind)?;
        Ok(Self::new(kind, len))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!("1M".parse::<Op>().unwrap(), Op::new(Kind::Match, 1));
        assert_eq!("13N".parse::<Op>().unwrap(), Op::new(Kind::Skip, 13));
        assert_eq!("144S".parse::<Op>().unwrap(), Op::new(Kind::SoftClip, 144));
    }
}
