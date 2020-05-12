mod key;

pub use self::key::Key;

use std::{error, fmt, ops::Deref, str::FromStr};

const DELIMITER: char = ':';

#[derive(Debug, Default)]
pub struct Format(Vec<Key>);

impl Deref for Format {
    type Target = [Key];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid format: {}", self.0)
    }
}

impl FromStr for Format {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.split(DELIMITER)
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map(Format)
            .map_err(|_| ParseError(s.into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Format = "GT".parse()?;
        assert_eq!(actual.len(), 1);

        let actual: Format = "GT:GQ:DP:HQ".parse()?;
        assert_eq!(actual.len(), 4);

        Ok(())
    }
}
