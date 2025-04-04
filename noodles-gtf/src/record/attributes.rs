mod field;

use std::{borrow::Cow, fmt, io};

use bstr::{BStr, BString, ByteSlice};
use indexmap::IndexMap;
use noodles_gff as gff;

use self::field::{parse_field, Value};

/// GTF record attributes.
pub struct Attributes<'r>(IndexMap<&'r BStr, Value<'r>>);

impl<'r> Attributes<'r> {
    pub(super) fn try_new(src: &'r [u8]) -> io::Result<Self> {
        parse_attributes(src).map(Self)
    }

    /// Returns whether there are any attributes.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the value of the given key.
    pub fn get(&self, key: &[u8]) -> Option<io::Result<&Value<'r>>> {
        self.0.get(key.as_bstr()).map(Ok)
    }

    /// Returns an iterator over key-value pairs.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<(&'r BStr, &Value<'r>)>> + '_ {
        self.0.iter().map(|(k, v)| Ok((*k, v)))
    }
}

impl fmt::Debug for Attributes<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut formatter = f.debug_map();

        for result in self.iter() {
            let (tag, value) = result.map_err(|_| fmt::Error)?;
            formatter.entry(&tag, value);
        }

        formatter.finish()
    }
}

impl gff::feature::record::Attributes for Attributes<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get(
        &self,
        tag: &[u8],
    ) -> Option<io::Result<gff::feature::record::attributes::field::Value<'_>>> {
        self.get(tag).map(|result| result.map(|value| value.into()))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<
                Item = io::Result<(
                    Cow<'_, BStr>,
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

fn parse_attributes(mut src: &[u8]) -> io::Result<IndexMap<&BStr, Value<'_>>> {
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

const BACKSLASH: u8 = b'\\';

fn escape_decode(s: &[u8]) -> io::Result<Cow<'_, BStr>> {
    if s.contains(&BACKSLASH) {
        unescape_string(s).map(Cow::from)
    } else {
        Ok(Cow::from(s.as_bstr()))
    }
}

fn unescape_string(s: &[u8]) -> io::Result<BString> {
    const QUOTATION_MARK: u8 = b'"';

    enum State {
        Normal,
        Escape,
    }

    let mut dst = Vec::with_capacity(s.len());
    let mut state = State::Normal;

    for c in s.iter().copied() {
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

    Ok(dst.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() -> io::Result<()> {
        assert!(Attributes::try_new(b"")?.is_empty());
        assert!(!Attributes::try_new(b"id 0;")?.is_empty());
        assert!(!Attributes::try_new(br#"id 0; name "ndls";"#)?.is_empty());
        Ok(())
    }

    #[test]
    fn test_get() -> io::Result<()> {
        let attributes = Attributes::try_new(br#"id 0; name "ndls";"#)?;

        assert_eq!(
            attributes.get(b"id").transpose()?,
            Some(&Value::String(Cow::from(BStr::new("0"))))
        );
        assert_eq!(
            attributes.get(b"name").transpose()?,
            Some(&Value::String(Cow::from(BStr::new("ndls"))))
        );

        assert!(attributes.get(b"comment").transpose()?.is_none());

        Ok(())
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        fn f(src: &[u8], expected: &[(&BStr, &Value<'_>)]) -> io::Result<()> {
            let attributes = Attributes::try_new(src)?;
            let actual = attributes.iter().collect::<io::Result<Vec<_>>>()?;
            assert_eq!(actual, expected);
            Ok(())
        }

        f(b"", &[])?;
        f(
            b"id 0;",
            &[(BStr::new("id"), &Value::String(Cow::from(BStr::new("0"))))],
        )?;
        f(
            b"id 0",
            &[(BStr::new("id"), &Value::String(Cow::from(BStr::new("0"))))],
        )?;
        f(
            br#"id 0; name "ndls";"#,
            &[
                (BStr::new("id"), &Value::String(Cow::from(BStr::new("0")))),
                (
                    BStr::new("name"),
                    &Value::String(Cow::from(BStr::new("ndls"))),
                ),
            ],
        )?;
        f(
            br#"id 0;name "ndls";"#,
            &[
                (BStr::new("id"), &Value::String(Cow::from(BStr::new("0")))),
                (
                    BStr::new("name"),
                    &Value::String(Cow::from(BStr::new("ndls"))),
                ),
            ],
        )?;
        f(
            br#"id 0;  name "ndls";  "#,
            &[
                (BStr::new("id"), &Value::String(Cow::from(BStr::new("0")))),
                (
                    BStr::new("name"),
                    &Value::String(Cow::from(BStr::new("ndls"))),
                ),
            ],
        )?;

        Ok(())
    }

    #[test]
    fn test_escape_decode() -> io::Result<()> {
        assert_eq!(escape_decode(b"")?, Cow::from(BStr::new("")));
        assert_eq!(escape_decode(b"ndls")?, Cow::from(BStr::new("ndls")));
        assert_eq!(escape_decode(br"nd\\ls")?, Cow::from(BStr::new(r"nd\ls")));
        assert_eq!(
            escape_decode(br#"nd\"ls\""#)?,
            Cow::from(BStr::new(r#"nd"ls""#))
        );

        assert!(matches!(
            escape_decode(br#"nd\ls"#),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
