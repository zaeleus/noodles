use std::{error, fmt, str::FromStr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    Info,
}

impl AsRef<str> for Kind {
    fn as_ref(&self) -> &str {
        match self {
            Self::Info => "INFO",
        }
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid record kind: expected {{INFO}}, got {}", self.0)
    }
}

impl FromStr for Kind {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "INFO" => Ok(Self::Info),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("INFO".parse::<Kind>()?, Kind::Info);

        assert!("".parse::<Kind>().is_err());
        assert!("##INFO".parse::<Kind>().is_err());

        Ok(())
    }
}
