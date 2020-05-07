use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    FileFormat,
    Info,
    Filter,
    Format,
    AlternativeAllele,
}

impl AsRef<str> for Kind {
    fn as_ref(&self) -> &str {
        match self {
            Self::FileFormat => "fileformat",
            Self::Info => "INFO",
            Self::Filter => "FILTER",
            Self::Format => "FORMAT",
            Self::AlternativeAllele => "ALT",
        }
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid record kind: expected {{fileformat, INFO, FILTER, FORMAT}}, got {}",
            self.0
        )
    }
}

impl FromStr for Kind {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "fileformat" => Ok(Self::FileFormat),
            "INFO" => Ok(Self::Info),
            "FILTER" => Ok(Self::Filter),
            "FORMAT" => Ok(Self::Format),
            "ALT" => Ok(Self::AlternativeAllele),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("fileformat".parse::<Kind>()?, Kind::FileFormat);
        assert_eq!("INFO".parse::<Kind>()?, Kind::Info);
        assert_eq!("FILTER".parse::<Kind>()?, Kind::Filter);
        assert_eq!("FORMAT".parse::<Kind>()?, Kind::Format);
        assert_eq!("ALT".parse::<Kind>()?, Kind::AlternativeAllele);

        assert!("".parse::<Kind>().is_err());
        assert!("##INFO".parse::<Kind>().is_err());

        Ok(())
    }
}
