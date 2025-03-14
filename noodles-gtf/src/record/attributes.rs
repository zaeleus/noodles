mod field;

use std::io;

use indexmap::IndexMap;

use self::field::{parse_field, Value};

/// GTF record attributes.
pub struct Attributes<'r>(IndexMap<&'r str, Value<'r>>);

impl<'r> Attributes<'r> {
    pub(super) fn try_new(src: &'r str) -> io::Result<Self> {
        parse_attributes(src).map(Self)
    }

    /// Returns whether there are any attributes.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the value of the given key.
    pub fn get(&self, key: &str) -> Option<io::Result<&Value<'r>>> {
        self.0.get(key).map(Ok)
    }

    /// Returns an iterator over key-value pairs.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<(&'r str, &Value<'r>)>> + '_ {
        self.0.iter().map(|(k, v)| Ok((*k, v)))
    }
}

fn parse_attributes(mut src: &str) -> io::Result<IndexMap<&str, Value<'_>>> {
    use indexmap::map::Entry;

    let mut map = IndexMap::new();

    while !src.is_empty() {
        let (key, raw_value) = parse_field(&mut src)?;

        match map.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(Value::String(raw_value));
            }
            Entry::Occupied(mut entry) => entry.get_mut().push(raw_value),
        }
    }

    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() -> io::Result<()> {
        assert!(Attributes::try_new("")?.is_empty());
        assert!(!Attributes::try_new("id 0;")?.is_empty());
        assert!(!Attributes::try_new(r#"id 0; name "ndls";"#)?.is_empty());
        Ok(())
    }

    #[test]
    fn test_get() -> io::Result<()> {
        let attributes = Attributes::try_new(r#"id 0; name "ndls";"#)?;

        assert_eq!(attributes.get("id").transpose()?, Some(&Value::String("0")));
        assert_eq!(
            attributes.get("name").transpose()?,
            Some(&Value::String("ndls"))
        );

        assert!(attributes.get("comment").transpose()?.is_none());

        Ok(())
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        fn f(src: &str, expected: &[(&str, &Value<'_>)]) -> io::Result<()> {
            let attributes = Attributes::try_new(src)?;
            let actual = attributes.iter().collect::<io::Result<Vec<_>>>()?;
            assert_eq!(actual, expected);
            Ok(())
        }

        f("", &[])?;
        f("id 0;", &[("id", &Value::String("0"))])?;
        f("id 0", &[("id", &Value::String("0"))])?;
        f(
            r#"id 0; name "ndls";"#,
            &[
                ("id", &Value::String("0")),
                ("name", &Value::String("ndls")),
            ],
        )?;
        f(
            r#"id 0;name "ndls";"#,
            &[
                ("id", &Value::String("0")),
                ("name", &Value::String("ndls")),
            ],
        )?;
        f(
            r#"id 0;  name "ndls";  "#,
            &[
                ("id", &Value::String("0")),
                ("name", &Value::String("ndls")),
            ],
        )?;

        Ok(())
    }
}
