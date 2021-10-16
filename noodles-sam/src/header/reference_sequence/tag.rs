//! SAM header reference sequence tag.

use std::{
    convert::TryFrom,
    error,
    fmt::{self, Write},
    str::FromStr,
};

const LENGTH: usize = 2;

/// A SAM header reference sequence tag.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
    /// Reference sequence name (`SN`).
    Name,
    /// Reference sequence length (`LN`).
    Length,
    /// Alternate locus (`AH`).
    AlternativeLocus,
    /// Alternate reference sequence names (`AN`).
    AlternativeNames,
    /// Genome assembly ID (`AS`).
    AssemblyId,
    /// Description (`DS`).
    Description,
    /// MD5 checksum of the reference sequence (`M5`).
    Md5Checksum,
    /// Species (`SP`).
    Species,
    /// Molecule topology (`TP`).
    MoleculeTopology,
    /// URI of the reference sequence (`UR`).
    Uri,
    /// Any other reference sequence tag.
    Other([u8; 2]),
}

impl AsRef<[u8; LENGTH]> for Tag {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Self::Name => b"SN",
            Self::Length => b"LN",
            Self::AlternativeLocus => b"AH",
            Self::AlternativeNames => b"AN",
            Self::AssemblyId => b"AS",
            Self::Description => b"DS",
            Self::Md5Checksum => b"M5",
            Self::Species => b"SP",
            Self::MoleculeTopology => b"TP",
            Self::Uri => b"UR",
            Self::Other(tag) => tag,
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

/// An error returned when a raw SAM header reference sequence tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Tag {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bytes = s.as_bytes();

        match bytes.len() {
            0 => Err(ParseError::Empty),
            2 => Self::try_from([bytes[0], bytes[1]]),
            _ => Err(ParseError::Invalid),
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Tag {
    type Error = ParseError;

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match &b {
            b"SN" => Ok(Self::Name),
            b"LN" => Ok(Self::Length),
            b"AH" => Ok(Self::AlternativeLocus),
            b"AN" => Ok(Self::AlternativeNames),
            b"AS" => Ok(Self::AssemblyId),
            b"DS" => Ok(Self::Description),
            b"M5" => Ok(Self::Md5Checksum),
            b"SP" => Ok(Self::Species),
            b"TP" => Ok(Self::MoleculeTopology),
            b"UR" => Ok(Self::Uri),
            _ => Ok(Self::Other(b)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Tag::Name.to_string(), "SN");
        assert_eq!(Tag::Length.to_string(), "LN");
        assert_eq!(Tag::AlternativeLocus.to_string(), "AH");
        assert_eq!(Tag::AlternativeNames.to_string(), "AN");
        assert_eq!(Tag::AssemblyId.to_string(), "AS");
        assert_eq!(Tag::Description.to_string(), "DS");
        assert_eq!(Tag::Md5Checksum.to_string(), "M5");
        assert_eq!(Tag::Species.to_string(), "SP");
        assert_eq!(Tag::MoleculeTopology.to_string(), "TP");
        assert_eq!(Tag::Uri.to_string(), "UR");
        assert_eq!(Tag::Other([b'N', b'D']).to_string(), "ND");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("SN".parse(), Ok(Tag::Name));
        assert_eq!("LN".parse(), Ok(Tag::Length));
        assert_eq!("AH".parse(), Ok(Tag::AlternativeLocus));
        assert_eq!("AN".parse(), Ok(Tag::AlternativeNames));
        assert_eq!("AS".parse(), Ok(Tag::AssemblyId));
        assert_eq!("DS".parse(), Ok(Tag::Description));
        assert_eq!("M5".parse(), Ok(Tag::Md5Checksum));
        assert_eq!("SP".parse(), Ok(Tag::Species));
        assert_eq!("TP".parse(), Ok(Tag::MoleculeTopology));
        assert_eq!("UR".parse(), Ok(Tag::Uri));
        assert_eq!("ND".parse(), Ok(Tag::Other([b'N', b'D'])));

        assert_eq!("".parse::<Tag>(), Err(ParseError::Empty));
        assert_eq!("NDL".parse::<Tag>(), Err(ParseError::Invalid));
    }
}
