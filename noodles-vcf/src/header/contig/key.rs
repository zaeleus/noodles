use std::{error, fmt, str::FromStr};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Key {
    Id,
    Other(String),
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid contig key: {}", self.0)
    }
}

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            "ID" => Ok(Self::Id),
            _ => Ok(Self::Other(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("ID".parse::<Key>()?, Key::Id);
        assert_eq!("length".parse::<Key>()?, Key::Other(String::from("length")));
        assert_eq!(
            "assembly".parse::<Key>()?,
            Key::Other(String::from("assembly"))
        );
        assert_eq!("md5".parse::<Key>()?, Key::Other(String::from("md5")));
        assert_eq!(
            "species".parse::<Key>()?,
            Key::Other(String::from("species"))
        );
        assert_eq!(
            "taxonomy".parse::<Key>()?,
            Key::Other(String::from("taxonomy"))
        );
        assert_eq!("URL".parse::<Key>()?, Key::Other(String::from("URL")));

        assert!("".parse::<Key>().is_err());

        Ok(())
    }
}
