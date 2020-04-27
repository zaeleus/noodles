use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    Header,
    ReferenceSequence,
    ReadGroup,
    Program,
    Comment,
}

impl AsRef<str> for Kind {
    fn as_ref(&self) -> &str {
        match self {
            Self::Header => "HD",
            Self::ReferenceSequence => "SQ",
            Self::ReadGroup => "RG",
            Self::Program => "PG",
            Self::Comment => "CO",
        }
    }
}

impl fmt::Display for Kind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "@{}", self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid record kind: expected @{{HD, SQ, RG, PG, CO}}, got {}",
            self.0
        )
    }
}

impl FromStr for Kind {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "@HD" => Ok(Self::Header),
            "@SQ" => Ok(Self::ReferenceSequence),
            "@RG" => Ok(Self::ReadGroup),
            "@PG" => Ok(Self::Program),
            "@CO" => Ok(Self::Comment),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("@HD".parse::<Kind>()?, Kind::Header);
        assert_eq!("@SQ".parse::<Kind>()?, Kind::ReferenceSequence);
        assert_eq!("@RG".parse::<Kind>()?, Kind::ReadGroup);
        assert_eq!("@PG".parse::<Kind>()?, Kind::Program);
        assert_eq!("@CO".parse::<Kind>()?, Kind::Comment);

        assert!("@NO".parse::<Kind>().is_err());
        assert!("HD".parse::<Kind>().is_err());

        Ok(())
    }
}
