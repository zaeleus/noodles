use std::{error, fmt, str::FromStr};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Key {
    Id,
    Length,
    Other(String),
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => "ID",
            Self::Length => "length",
            Self::Other(s) => s,
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
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
            "length" => Ok(Self::Length),
            _ => Ok(Self::Other(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Key::Id.to_string(), "ID");
        assert_eq!(Key::Length.to_string(), "length");
        assert_eq!(Key::Other(String::from("md5")).to_string(), "md5");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("ID".parse::<Key>()?, Key::Id);
        assert_eq!("length".parse::<Key>()?, Key::Length);
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
