use std::{error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Id(Option<String>);

impl AsRef<str> for Id {
    fn as_ref(&self) -> &str {
        self.0.as_deref().unwrap_or(MISSING_FIELD)
    }
}

impl Deref for Id {
    type Target = Option<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Id {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid id: {}", self.0)
    }
}

impl FromStr for Id {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            MISSING_FIELD => Ok(Self(None)),
            _ => Ok(Self(Some(s.into()))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Id(Some(String::from("r0"))).to_string(), "r0");
        assert_eq!(Id(None).to_string(), ".");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert!(".".parse::<Id>()?.is_none());
        assert!("rs13".parse::<Id>()?.is_some());
        assert!("".parse::<Id>().is_err());
        Ok(())
    }
}
