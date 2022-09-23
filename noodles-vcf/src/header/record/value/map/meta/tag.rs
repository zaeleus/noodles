use std::str::FromStr;

use crate::header::record::value::map::tag::{self, ID, NUMBER, TYPE};

const VALUES: &str = "Values";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Standard {
    Id,
    Type,
    Number,
    Values,
}

impl tag::Standard for Standard {}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => tag::ID,
            Self::Type => tag::TYPE,
            Self::Number => tag::NUMBER,
            Self::Values => VALUES,
        }
    }
}

impl FromStr for Standard {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            ID => Ok(Self::Id),
            TYPE => Ok(Self::Type),
            NUMBER => Ok(Self::Number),
            VALUES => Ok(Self::Values),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_as_ref_str_for_standard() {
        assert_eq!(Standard::Id.as_ref(), "ID");
        assert_eq!(Standard::Type.as_ref(), "Type");
        assert_eq!(Standard::Number.as_ref(), "Number");
        assert_eq!(Standard::Values.as_ref(), "Values");
    }
}
