pub mod field;

pub use self::field::Field;

use std::{error, fmt, ops::Deref, str::FromStr};

const DELIMITER: char = ';';

#[derive(Debug, Default)]
pub struct Info(Vec<Field>);

impl Deref for Info {
    type Target = [Field];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid info: {}", self.0)
    }
}

impl FromStr for Info {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.split(DELIMITER)
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map(Info)
            .map_err(|_| ParseError(s.into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Info = "NS=2".parse()?;
        assert_eq!(actual.len(), 1);

        let actual: Info = "NS=2;AF=0.333,0.667".parse()?;
        assert_eq!(actual.len(), 2);

        Ok(())
    }
}
