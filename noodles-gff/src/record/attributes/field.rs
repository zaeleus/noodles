//! GFF record attributes field.

mod tag;
mod value;

use std::{borrow::Cow, io};

use bstr::{BStr, ByteSlice};

pub use self::value::Value;
use self::{tag::parse_tag, value::parse_value};

pub(super) fn parse_field<'a>(src: &mut &'a [u8]) -> io::Result<(Cow<'a, BStr>, Value<'a>)> {
    const DELIMITER: u8 = b';';
    const SEPARATOR: u8 = b'=';

    let (buf, rest) = split_once(src, DELIMITER).unwrap_or_else(|| src.split_at(src.len()));

    *src = rest;

    if let Some((t, v)) = split_once(buf, SEPARATOR) {
        let tag = parse_tag(t);
        let value = parse_value(v);
        Ok((tag, value))
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
        let mut src = &b"ID=1;Name=ndls"[..];
        assert_eq!(
            parse_field(&mut src)?,
            (
                Cow::from(BStr::new("ID")),
                Value::String(Cow::from(BStr::new("1")))
            )
        );
        assert_eq!(
            parse_field(&mut src)?,
            (
                Cow::from(BStr::new("Name")),
                Value::String(Cow::from(BStr::new("ndls")))
            )
        );
        assert!(src.is_empty());

        let mut src = &b"ID"[..];
        assert!(matches!(
            parse_field(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
