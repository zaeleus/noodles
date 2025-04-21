//! GFF record attributes.

pub mod field;

use std::{borrow::Cow, fmt, io, iter};

use bstr::BStr;

use self::field::{next_field, parse_field, parse_tag, parse_value, split_field, Value};

/// GFF record attributes.
pub struct Attributes<'a>(&'a [u8]);

impl<'a> Attributes<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }

    /// Returns whether there are any attributes.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the value of the given tag.
    pub fn get(&self, tag: &[u8]) -> Option<io::Result<Value<'_>>> {
        let mut src = self.0;

        while let Some(s) = next_field(&mut src) {
            let (t, v) = match split_field(s) {
                Ok(srcs) => srcs,
                Err(e) => return Some(Err(e)),
            };

            if parse_tag(t).as_ref() == tag {
                return Some(Ok(parse_value(v)));
            }
        }

        None
    }

    /// Returns an iterator over all tag-value pairs.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<(Cow<'_, BStr>, Value<'_>)>> {
        let mut src = self.0;
        iter::from_fn(move || next_field(&mut src).map(parse_field))
    }
}

impl AsRef<[u8]> for Attributes<'_> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl fmt::Debug for Attributes<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut formatter = f.debug_map();

        for result in self.iter() {
            let (tag, value) = result.map_err(|_| fmt::Error)?;
            formatter.entry(&tag, &value);
        }

        formatter.finish()
    }
}

impl crate::feature::record::Attributes for Attributes<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get(
        &self,
        tag: &[u8],
    ) -> Option<io::Result<crate::feature::record::attributes::field::Value<'_>>> {
        self.get(tag).map(|result| result.map(|value| value.into()))
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<
                Item = io::Result<(
                    Cow<'_, BStr>,
                    crate::feature::record::attributes::field::Value<'_>,
                )>,
            > + '_,
    > {
        Box::new(
            self.iter()
                .map(|result| result.map(|(tag, value)| (tag, value.into()))),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        let attributes = Attributes::new(b"");
        assert!(attributes.is_empty());

        let attributes = Attributes::new(b"ID=1;Name=ndls");
        assert!(!attributes.is_empty());
    }

    #[test]
    fn test_get() -> io::Result<()> {
        let attributes = Attributes::new(b"ID=1;Name%3F=ndls");

        assert_eq!(
            attributes.get(b"ID").transpose()?,
            Some(Value::String(Cow::from(BStr::new("1"))))
        );

        assert_eq!(
            attributes.get(b"Name?").transpose()?,
            Some(Value::String(Cow::from(BStr::new("ndls"))))
        );

        assert!(attributes.get(b"comment").is_none());

        Ok(())
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        let attributes = Attributes::new(b"");
        assert!(attributes.iter().next().is_none());

        let attributes = Attributes::new(b"ID=1;Name=ndls");
        let actual: Vec<_> = attributes.iter().collect::<Result<_, _>>()?;
        let expected = vec![
            (
                Cow::from(BStr::new("ID")),
                Value::String(Cow::from(BStr::new("1"))),
            ),
            (
                Cow::from(BStr::new("Name")),
                Value::String(Cow::from(BStr::new("ndls"))),
            ),
        ];
        assert_eq!(actual, expected);

        Ok(())
    }
}
