//! GTF record attributes field.

mod value;

use std::io;

use bstr::{BStr, ByteSlice};

pub use self::value::Value;

const SEPARATOR: u8 = b' ';
const DOUBLE_QUOTES: u8 = b'"';
const TERMINATOR: u8 = b';';

pub(super) fn parse_field<'a>(src: &mut &'a [u8]) -> io::Result<(&'a BStr, &'a BStr)> {
    let key = parse_key(src)?;
    let value = parse_value(src)?;
    maybe_consume_terminator(src);
    Ok((key, value))
}

fn parse_key<'a>(src: &mut &'a [u8]) -> io::Result<&'a BStr> {
    let Some(i) = src.iter().position(|c| *c == SEPARATOR) else {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    };

    let (buf, rest) = src.split_at(i);
    *src = &rest[1..];
    Ok(buf.as_bstr())
}

fn parse_value<'a>(src: &mut &'a [u8]) -> io::Result<&'a BStr> {
    if let Some(rest) = src.strip_prefix(&[DOUBLE_QUOTES]) {
        *src = rest;
        parse_string(src)
    } else {
        parse_raw_value(src)
    }
}

fn parse_string<'a>(src: &mut &'a [u8]) -> io::Result<&'a BStr> {
    let Some(i) = src.iter().position(|c| *c == DOUBLE_QUOTES) else {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    };

    let (buf, rest) = src.split_at(i);
    *src = &rest[1..];
    Ok(buf.as_bstr())
}

fn parse_raw_value<'a>(src: &mut &'a [u8]) -> io::Result<&'a BStr> {
    let i = src
        .iter()
        .position(|c| *c == TERMINATOR)
        .unwrap_or(src.len());

    let (buf, rest) = src.split_at(i);
    *src = rest;
    Ok(buf.as_bstr())
}

fn maybe_consume_terminator(src: &mut &[u8]) {
    *src = src.trim_ascii_start();

    if let Some(rest) = src.strip_prefix(&[TERMINATOR]) {
        *src = rest.trim_ascii_start();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() -> io::Result<()> {
        fn t(mut src: &[u8], (expected_key, expected_value): (&[u8], &[u8])) -> io::Result<()> {
            assert_eq!(
                parse_field(&mut src)?,
                (BStr::new(expected_key), BStr::new(expected_value))
            );

            Ok(())
        }

        t(br#"id "0""#, (b"id", b"0"))?;
        t(b"id 0", (b"id", b"0"))?;
        t(br#"id "0";"#, (b"id", b"0"))?;
        t(b"id 0;", (b"id", b"0"))?;
        t(br#"id "0" ; "#, (b"id", b"0"))?;
        t(br#"id "0;1";"#, (b"id", b"0;1"))?;

        assert!(matches!(
            parse_field(&mut &b""[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            parse_field(&mut &b"id"[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
