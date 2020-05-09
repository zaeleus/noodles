use std::{error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

const DELIMITER: char = ',';

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlternateBases(Vec<String>);

impl Deref for AlternateBases {
    type Target = [String];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid alternate bases: {}", self.0)
    }
}

impl FromStr for AlternateBases {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            MISSING_FIELD => Ok(AlternateBases::default()),
            _ => Ok(AlternateBases(
                s.split(DELIMITER).map(|s| s.into()).collect(),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert!(".".parse::<AlternateBases>()?.is_empty());

        let expected = AlternateBases(vec![String::from("G")]);
        assert_eq!("G".parse::<AlternateBases>()?, expected);

        let expected = AlternateBases(vec![String::from("G"), String::from("T")]);
        assert_eq!("G,T".parse::<AlternateBases>()?, expected);

        let expected = AlternateBases(vec![String::from("<DUP>")]);
        assert_eq!("<DUP>".parse::<AlternateBases>()?, expected);

        let expected = AlternateBases(vec![String::from("]sq0:5]A")]);
        assert_eq!("]sq0:5]A".parse::<AlternateBases>()?, expected);

        let expected = AlternateBases(vec![String::from("C[sq1:13[")]);
        assert_eq!("C[sq1:13[".parse::<AlternateBases>()?, expected);

        assert!("".parse::<AlternateBases>().is_err());

        Ok(())
    }
}
