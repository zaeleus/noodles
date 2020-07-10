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

impl fmt::Display for FilterStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing => f.write_str(MISSING_FIELD),
            Self::Pass => f.write_str(PASS_STATUS),
            Self::Fail(ids) => {
                for (i, id) in ids.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{}", DELIMITER)?;
                    }

                    f.write_str(id)?;
                }

                Ok(())
            }
        }
    }
}

/// An error returned when a raw VCF filter status fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
        }
    }
}

impl FromStr for FilterStatus {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
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
    fn test_fmt() {
        assert_eq!(FilterStatus::Missing.to_string(), ".");
        assert_eq!(FilterStatus::Pass.to_string(), "PASS");

        let status = FilterStatus::Fail(vec![String::from("q10")]);
        assert_eq!(status.to_string(), "q10");

        let status = FilterStatus::Fail(vec![String::from("q10"), String::from("s50")]);
        assert_eq!(status.to_string(), "q10,s50");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(".".parse(), Ok(FilterStatus::Missing));
        assert_eq!("PASS".parse(), Ok(FilterStatus::Pass));

        assert_eq!(
            "q10".parse(),
            Ok(FilterStatus::Fail(vec![String::from("q10")]))
        );

        assert_eq!(
            "q10,s50".parse(),
            Ok(FilterStatus::Fail(vec![
                String::from("q10"),
                String::from("s50")
            ]))
        );

        assert_eq!("".parse::<FilterStatus>(), Err(ParseError::Empty));
    }
}
