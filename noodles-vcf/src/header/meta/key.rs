//! VCF header meta record key.

use std::fmt;

/// A VCF header meta record key.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Key {
    /// (`ID`).
    Id,
    /// (`Type`).
    Type,
    /// (`Number`).
    Number,
    /// (`Values`).
    Values,
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => "ID",
            Self::Type => "Type",
            Self::Number => "Number",
            Self::Values => "Values",
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Key::Id.to_string(), "ID");
        assert_eq!(Key::Type.to_string(), "Type");
        assert_eq!(Key::Number.to_string(), "Number");
        assert_eq!(Key::Values.to_string(), "Values");
    }
}
