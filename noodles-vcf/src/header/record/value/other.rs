use std::fmt;

use super::{map, Map};

/// A VCF header record value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Other {
    /// A string.
    String(String),
    /// A map.
    Map(String, Map<map::Other>),
}

impl fmt::Display for Other {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::String(s) => s.fmt(f),
            Self::Map(id, map) => write!(f, "<ID={id}{map}>"),
        }
    }
}

impl From<&str> for Other {
    fn from(s: &str) -> Self {
        Self::String(s.into())
    }
}

impl From<String> for Other {
    fn from(s: String) -> Self {
        Self::String(s)
    }
}

impl From<(String, Map<map::Other>)> for Other {
    fn from((id, map): (String, Map<map::Other>)) -> Self {
        Self::Map(id, map)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Other::from("VCFv4.3").to_string(), "VCFv4.3");
        assert_eq!(
            Other::Map(String::from("0"), Map::<map::Other>::new()).to_string(),
            "<ID=0>"
        );
    }

    #[test]
    fn test_from_str_for_other() {
        let s = "noodles";
        assert_eq!(Other::from(s), Other::String(s.into()));
    }

    #[test]
    fn test_from_string_for_other() {
        let s = String::from("noodles");
        assert_eq!(Other::from(s.clone()), Other::String(s));
    }

    #[test]
    fn test_from_map_other_for_value() {
        let id = String::from("0");
        let map = Map::<map::Other>::new();
        assert_eq!(Other::from((id.clone(), map.clone())), Other::Map(id, map));
    }
}
