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

    /// Returns the value of the given key.
    pub fn get(&self, key: &str) -> Option<io::Result<&str>> {
        for result in self.iter() {
            match result {
                Ok((k, value)) => {
                    if k == key {
                        return Some(Ok(value));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
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
    fn test_get() -> io::Result<()> {
        let attributes = Attributes::new(r#"id 0; name "ndls";"#);

        assert_eq!(attributes.get("id").transpose()?, Some("0"));
        assert_eq!(attributes.get("name").transpose()?, Some("ndls"));

        assert!(attributes.get("comment").transpose()?.is_none());

        Ok(())
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        fn f(src: &str, expected: &[(&str, &str)]) -> io::Result<()> {
            let actual = Attributes::new(src)
                .iter()
                .collect::<io::Result<Vec<_>>>()?;

            assert_eq!(actual, expected);

            Ok(())
        }

        f("", &[])?;
        f("id 0;", &[("id", "0")])?;
        f("id 0", &[("id", "0")])?;
        f(r#"id 0; name "ndls";"#, &[("id", "0"), ("name", "ndls")])?;
        f(r#"id 0;name "ndls";"#, &[("id", "0"), ("name", "ndls")])?;
        f(r#"id 0;  name "ndls";  "#, &[("id", "0"), ("name", "ndls")])?;

        Ok(())
    }
}
