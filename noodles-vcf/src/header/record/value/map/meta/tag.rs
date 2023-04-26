use std::str::FromStr;

use crate::header::record::value::map;

pub(super) type StandardTag = Standard;

/// A VCF header meta map tag.
pub type Tag = map::tag::Tag<StandardTag>;

// For some reason, using the `Tag` type alias produces a `nontrivial_structural_match` warning
// when pattern matching, so it's avoided here.
pub(crate) const ID: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Id);
pub(super) const TYPE: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Type);
pub(super) const NUMBER: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Number);
pub(super) const VALUES: Tag = map::tag::Tag::<StandardTag>::Standard(StandardTag::Values);

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
            Self::Id => "ID",
            Self::Type => "Type",
            Self::Number => "Number",
            Self::Values => "Values",
        }
    }
}

impl FromStr for Standard {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "ID" => Ok(Self::Id),
            "Type" => Ok(Self::Type),
            "Number" => Ok(Self::Number),
            "Values" => Ok(Self::Values),
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
