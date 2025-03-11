mod field;

use std::{io, iter};

use self::field::parse_field;

/// GTF record attributes.
pub struct Attributes<'r>(&'r str);

impl<'r> Attributes<'r> {
    pub(super) fn new(src: &'r str) -> Self {
        Self(src)
    }

    /// Returns whether there are any attributes.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns an iterator over key-value pairs.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<(&'r str, &'r str)>> {
        let mut src = self.0;

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(parse_field(&mut src))
            }
        })
    }
}

impl AsRef<str> for Attributes<'_> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        assert!(Attributes::new("").is_empty());
        assert!(!Attributes::new("id 0;").is_empty());
        assert!(!Attributes::new(r#"id 0; name "ndls";"#).is_empty());
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        assert_eq!(
            Attributes::new("").iter().collect::<io::Result<Vec<_>>>()?,
            []
        );

        assert_eq!(
            Attributes::new("id 0;")
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [("id", "0")]
        );

        assert_eq!(
            Attributes::new(r#"id 0; name "ndls";"#)
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [("id", "0"), ("name", "ndls")]
        );

        Ok(())
    }
}
