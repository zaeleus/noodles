use std::{error, fmt, str::FromStr};

/// A SAM header header group order (`GO`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GroupOrder {
    /// Alignments are not grouped (`none`).
    None,
    /// Alignments are grouped by read name (`query`).
    Query,
    /// Alignments are grouped by reference sequence and position (`reference`).
    Reference,
}

impl AsRef<str> for GroupOrder {
    fn as_ref(&self) -> &str {
        match self {
            Self::None => "none",
            Self::Query => "query",
            Self::Reference => "reference",
        }
    }
}

impl Default for GroupOrder {
    fn default() -> Self {
        Self::None
    }
}

impl fmt::Display for GroupOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid group order: expected {{none, query, reference}}, got {}",
            self.0
        )
    }
}

impl FromStr for GroupOrder {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "none" => Ok(Self::None),
            "query" => Ok(Self::Query),
            "reference" => Ok(Self::Reference),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(GroupOrder::default(), GroupOrder::None);
    }

    #[test]
    fn test_fmt() {
        assert_eq!(GroupOrder::None.to_string(), "none");
        assert_eq!(GroupOrder::Query.to_string(), "query");
        assert_eq!(GroupOrder::Reference.to_string(), "reference");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("none".parse::<GroupOrder>()?, GroupOrder::None);
        assert_eq!("query".parse::<GroupOrder>()?, GroupOrder::Query);
        assert_eq!("reference".parse::<GroupOrder>()?, GroupOrder::Reference);

        assert!("".parse::<GroupOrder>().is_err());
        assert!("noodles".parse::<GroupOrder>().is_err());
        assert!("Query".parse::<GroupOrder>().is_err());

        Ok(())
    }
}
