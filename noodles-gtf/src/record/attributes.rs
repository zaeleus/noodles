mod field;

use std::{borrow::Cow, io};

use indexmap::IndexMap;
use noodles_gff as gff;

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

impl gff::feature::record::Attributes for Attributes<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get(
        &self,
        tag: &str,
    ) -> Option<io::Result<gff::feature::record::attributes::field::Value<'_>>> {
        self.get(tag).map(|result| result.map(|value| value.into()))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<
                Item = io::Result<(
                    Cow<'_, str>,
                    gff::feature::record::attributes::field::Value<'_>,
                )>,
            > + '_,
    > {
        Box::new(
            self.iter()
                .map(|result| result.map(|(key, value)| (Cow::from(key), value.into()))),
        )
    }
}

fn parse_attributes(mut src: &str) -> io::Result<IndexMap<&str, Value<'_>>> {
    use indexmap::map::Entry;

    let mut map = IndexMap::new();

    while !src.is_empty() {
        let (key, raw_value) = parse_field(&mut src)?;
        let decoded_raw_value = escape_decode(raw_value)?;

        match map.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(Value::String(decoded_raw_value));
            }
            Entry::Occupied(mut entry) => entry.get_mut().push(decoded_raw_value),
        }
    }

    Ok(map)
}

fn escape_decode(s: &str) -> io::Result<Cow<'_, str>> {
    const BACKSLASH: char = '\\';

    if s.contains(BACKSLASH) {
        unescape_string(s).map(Cow::from)
    } else {
        Ok(Cow::from(s))
    }
}

fn unescape_string(s: &str) -> io::Result<String> {
    const BACKSLASH: char = '\\';
    const QUOTATION_MARK: char = '"';

    enum State {
        Normal,
        Escape,
    }

    let mut dst = String::with_capacity(s.len());
    let mut state = State::Normal;

    for c in s.chars() {
        match state {
            State::Normal => {
                if c == BACKSLASH {
                    state = State::Escape;
                } else {
                    dst.push(c);
                }
            }
            State::Escape => {
                match c {
                    BACKSLASH | QUOTATION_MARK => dst.push(c),
                    _ => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "invalid escape sequence",
                        ))
                    }
                }

                state = State::Normal;
            }
        }
    }

    Ok(dst)
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

        assert_eq!(
            attributes.get("id").transpose()?,
            Some(&Value::String(Cow::from("0")))
        );
        assert_eq!(
            attributes.get("name").transpose()?,
            Some(&Value::String(Cow::from("ndls")))
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
        f("id 0;", &[("id", &Value::String(Cow::from("0")))])?;
        f("id 0", &[("id", &Value::String(Cow::from("0")))])?;
        f(
            r#"id 0; name "ndls";"#,
            &[
                ("id", &Value::String(Cow::from("0"))),
                ("name", &Value::String(Cow::from("ndls"))),
            ],
        )?;
        f(
            r#"id 0;name "ndls";"#,
            &[
                ("id", &Value::String(Cow::from("0"))),
                ("name", &Value::String(Cow::from("ndls"))),
            ],
        )?;
        f(
            r#"id 0;  name "ndls";  "#,
            &[
                ("id", &Value::String(Cow::from("0"))),
                ("name", &Value::String(Cow::from("ndls"))),
            ],
        )?;

        Ok(())
    }

    #[test]
    fn test_escape_decode() -> io::Result<()> {
        assert_eq!(escape_decode("")?, "");
        assert_eq!(escape_decode("ndls")?, "ndls");
        assert_eq!(escape_decode(r"nd\\ls")?, r"nd\ls");
        assert_eq!(escape_decode(r#"nd\"ls\""#)?, r#"nd"ls""#);

        assert!(matches!(
            escape_decode(r#"nd\ls"#),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
