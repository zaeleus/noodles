//! GFF record attributes field.

mod tag;
mod value;

use std::{borrow::Cow, io};

use bstr::{BStr, ByteSlice};

pub use self::value::Value;
pub(super) use self::{tag::parse_tag, value::parse_value};

pub(super) fn parse<'a>(t: &'a [u8], v: &'a [u8]) -> (Cow<'a, BStr>, Value<'a>) {
    (parse_tag(t), parse_value(v))
}

pub(super) fn next<'a>(src: &mut &'a [u8]) -> Option<io::Result<(&'a [u8], &'a [u8])>> {
    if src.is_empty() {
        return None;
    }

    let tag = match take_tag(src) {
        Ok(buf) => buf,
        Err(e) => return Some(Err(e)),
    };

    let value = take_value(src);

    Some(Ok((tag, value)))
}

fn take_tag<'a>(src: &mut &'a [u8]) -> io::Result<&'a [u8]> {
    const SEPARATOR: u8 = b'=';

    if let Some((buf, rest)) = split_once(src, SEPARATOR) {
        *src = rest;
        Ok(buf)
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidData, "invalid tag"))
    }
}

fn take_value<'a>(src: &mut &'a [u8]) -> &'a [u8] {
    const DELIMITER: u8 = b';';

    if let Some((buf, rest)) = split_once(src, DELIMITER) {
        *src = rest;
        buf
    } else {
        let (buf, rest) = src.split_at(src.len());
        *src = rest;
        buf
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
    fn test_next() -> io::Result<()> {
        let mut src = &b""[..];
        assert!(next(&mut src).is_none());

        let mut src = &b"ID=1;Name%3F=ndls"[..];
        assert_eq!(next(&mut src).transpose()?, Some((&b"ID"[..], &b"1"[..])));
        assert_eq!(
            next(&mut src).transpose()?,
            Some((&b"Name%3F"[..], &b"ndls"[..]))
        );

        let mut src = &b"ID"[..];
        assert!(matches!(
            next(&mut src),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
