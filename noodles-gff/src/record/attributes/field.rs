//! GFF record attributes field.

mod tag;
mod value;

use std::{borrow::Cow, io};

use bstr::{BStr, ByteSlice};

pub use self::value::Value;
use self::{tag::parse_tag, value::parse_value};

pub(super) fn parse_field(src: &[u8]) -> io::Result<(Cow<'_, BStr>, Value<'_>)> {
    split_field(src).map(|(t, v)| (parse_tag(t), parse_value(v)))
}

pub(super) fn next_field<'a>(src: &mut &'a [u8]) -> Option<&'a [u8]> {
    const DELIMITER: u8 = b';';

    let (buf, rest) = split_once(src, DELIMITER).unwrap_or_else(|| src.split_at(src.len()));

    *src = rest;

    if buf.is_empty() {
        None
    } else {
        Some(buf)
    }
}

fn split_field(src: &[u8]) -> io::Result<(&[u8], &[u8])> {
    const SEPARATOR: u8 = b'=';

    if let Some((t, v)) = split_once(src, SEPARATOR) {
        Ok((t, v))
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidData, "invalid field"))
    }
}

fn split_once(src: &[u8], n: u8) -> Option<(&[u8], &[u8])> {
    let i = src.iter().position(|b| *b == n)?;
    Some((&src[..i], &src[i + 1..]))
}

fn percent_decode(s: &[u8]) -> Cow<'_, BStr> {
    match Cow::from(percent_encoding::percent_decode(s)) {
        Cow::Borrowed(buf) => Cow::Borrowed(buf.as_bstr()),
        Cow::Owned(buf) => Cow::Owned(buf.into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() -> io::Result<()> {
        assert_eq!(
            parse_field(b"ID=1")?,
            (
                Cow::from(BStr::new("ID")),
                Value::String(Cow::from(BStr::new("1")))
            )
        );
        assert_eq!(
            parse_field(b"Name=ndls")?,
            (
                Cow::from(BStr::new("Name")),
                Value::String(Cow::from(BStr::new("ndls")))
            )
        );

        assert!(matches!(
            parse_field(b"ID"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
