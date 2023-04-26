use std::str::FromStr;

use crate::header::record::value::map::{
    self,
    tag::{ID, IDX},
};

pub(super) type StandardTag = Standard;

/// A VCF header contig map tag.
pub type Tag = map::tag::Tag<StandardTag>;

const LENGTH: &str = "length";
const MD5: &str = "md5";
const URL: &str = "URL";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Standard {
    Id,
    Length,
    Md5,
    Url,
    Idx,
}

impl map::tag::Standard for Standard {}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => ID,
            Self::Length => LENGTH,
            Self::Md5 => MD5,
            Self::Url => URL,
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
            MD5 => Ok(Self::Md5),
            URL => Ok(Self::Url),
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
        assert_eq!(Standard::Md5.as_ref(), "md5");
        assert_eq!(Standard::Url.as_ref(), "URL");
        assert_eq!(Standard::Idx.as_ref(), "IDX");
    }
}
