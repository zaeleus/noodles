//! BED record name.

use std::{fmt, ops::Deref, str::FromStr};

/// A BED record name.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Name(String);

impl Deref for Name {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self)
    }
}

impl FromStr for Name {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(s.into()))
    }
}

impl From<&str> for Name {
    fn from(s: &str) -> Self {
        Self::from(s.to_owned())
    }
}

impl From<String> for Name {
    fn from(s: String) -> Self {
        Self(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let name = Name(String::from("ndls1"));
        assert_eq!(name.to_string(), "ndls1");
    }

    #[test]
    fn test_from_str() {
        const MAX_LENGTH: usize = 255;

        assert_eq!("ndls1".parse(), Ok(Name(String::from("ndls1"))));
        assert_eq!(" ~".parse(), Ok(Name(String::from(" ~"))));

        assert_eq!("".parse::<Name>(), Ok(Name(String::new())));
        assert_eq!("🍜".parse::<Name>(), Ok(Name(String::from("🍜"))));

        let s = "n".repeat(MAX_LENGTH + 1);
        assert_eq!(s.parse::<Name>(), Ok(Name::from(s)));
    }
}
