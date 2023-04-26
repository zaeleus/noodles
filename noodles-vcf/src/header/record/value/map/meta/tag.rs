use std::str::FromStr;

use crate::header::record::value::map::{
    self,
    tag::{ID, NUMBER, TYPE},
};

pub(super) type StandardTag = Standard;

/// A VCF header meta map tag.
pub type Tag = map::tag::Tag<StandardTag>;

const VALUES: &str = "Values";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Standard {
    Id,
    Type,
    Number,
    Values,
}

impl map::tag::Standard for Standard {}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => map::tag::ID,
            Self::Type => map::tag::TYPE,
            Self::Number => map::tag::NUMBER,
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
