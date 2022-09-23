use std::str::FromStr;

use crate::header::record::value::map::tag::{self, ID, IDX};

const LENGTH: &str = "length";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Standard {
    Id,
    Length,
    Idx,
}

impl tag::Standard for Standard {}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => ID,
            Self::Length => LENGTH,
            Self::Idx => IDX,
        }
    }
}

impl FromStr for Standard {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            ID => Ok(Self::Id),
            LENGTH => Ok(Self::Length),
            IDX => Ok(Self::Idx),
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
        assert_eq!(Standard::Length.as_ref(), "length");
        assert_eq!(Standard::Idx.as_ref(), "IDX");
    }
}
