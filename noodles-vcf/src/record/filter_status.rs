use std::{error, fmt, str::FromStr};

use super::MISSING_FIELD;

const PASS_STATUS: &str = "PASS";
const DELIMITER: char = ',';

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum FilterStatus {
    Missing,
    Pass,
    Fail(Vec<String>),
}

impl Default for FilterStatus {
    fn default() -> Self {
        Self::Missing
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid filter status: {}", self.0)
    }
}

impl FromStr for FilterStatus {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError(s.into())),
            MISSING_FIELD => Ok(Self::Missing),
            PASS_STATUS => Ok(Self::Pass),
            _ => {
                let filters = s.split(DELIMITER).map(|t| t.into()).collect();
                Ok(Self::Fail(filters))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(FilterStatus::default(), FilterStatus::Missing);
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!(".".parse::<FilterStatus>()?, FilterStatus::Missing);
        assert_eq!("PASS".parse::<FilterStatus>()?, FilterStatus::Pass);

        let status = FilterStatus::Fail(vec![String::from("q10")]);
        assert_eq!("q10".parse::<FilterStatus>()?, status);

        let status = FilterStatus::Fail(vec![String::from("q10"), String::from("s50")]);
        assert_eq!("q10,s50".parse::<FilterStatus>()?, status);

        assert!("".parse::<FilterStatus>().is_err());

        Ok(())
    }
}
