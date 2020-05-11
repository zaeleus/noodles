use std::{error, fmt, str::FromStr};

const LEN: usize = 2;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Tag {
    Cigar,
    ReadGroup,
    AlignmentHitCount,
    Other(String),
}

impl AsRef<str> for Tag {
    fn as_ref(&self) -> &str {
        match self {
            Self::Cigar => "CG",
            Self::AlignmentHitCount => "NH",
            Self::ReadGroup => "RG",
            Self::Other(tag) => tag,
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid data tag: {}", self.0)
    }
}

impl FromStr for Tag {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "CG" => Ok(Tag::Cigar),
            "NH" => Ok(Tag::AlignmentHitCount),
            "RG" => Ok(Tag::ReadGroup),
            _ => {
                if s.len() == LEN {
                    Ok(Tag::Other(s.into()))
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
        assert_eq!(Tag::Cigar.to_string(), "CG");
        assert_eq!(Tag::AlignmentHitCount.to_string(), "NH");
        assert_eq!(Tag::ReadGroup.to_string(), "RG");
        assert_eq!(Tag::Other(String::from("ZN")).to_string(), "ZN");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("CG".parse::<Tag>()?, Tag::Cigar);
        assert_eq!("NH".parse::<Tag>()?, Tag::AlignmentHitCount);
        assert_eq!("RG".parse::<Tag>()?, Tag::ReadGroup);
        assert_eq!("ZN".parse::<Tag>()?, Tag::Other(String::from("ZN")));

        assert!("".parse::<Tag>().is_err());
        assert!("R".parse::<Tag>().is_err());
        assert!("RGP".parse::<Tag>().is_err());
        assert!("RGRP".parse::<Tag>().is_err());

        Ok(())
    }
}
