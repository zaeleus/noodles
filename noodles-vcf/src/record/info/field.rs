mod value;

use std::io;

use crate::{Header, header::record::value::map::info::Type, variant::record::info::field::Value};

const DELIMITER: u8 = b';';
const SEPARATOR: u8 = b'=';

pub(super) fn parse_value<'a>(
    header: &Header,
    key: &str,
    raw_value: Option<&'a str>,
) -> io::Result<Option<Value<'a>>> {
    use crate::header::record::value::map::info::definition::definition;

    const MISSING: &str = ".";

    let definition = header
        .infos()
        .get(key)
        .map(|info| (info.number(), info.ty()))
        .or_else(|| definition(header.file_format(), key).map(|(n, t, _)| (n, t)));

    if definition.is_none() && raw_value.is_none() {
        return Ok(Some(Value::Flag));
    }

    let (number, ty) = definition.unwrap_or_default();

    match raw_value {
        Some(MISSING) => Ok(None),
        Some(v) => value::parse_value(v, number, ty).map(Some),
        None if ty == Type::Flag => Ok(Some(Value::Flag)),
        None => Err(io::Error::new(io::ErrorKind::InvalidData, "missing value")),
    }
}

pub(super) fn next<'a>(src: &mut &'a str) -> Option<io::Result<(&'a str, Option<&'a str>)>> {
    if src.is_empty() {
        return None;
    }

    let (key, is_separated) = match read_key(src) {
        Ok((k, is_eof)) => (k, is_eof),
        Err(e) => return Some(Err(e)),
    };

    if !is_separated {
        return Some(Ok((key, None)));
    }

    let value = match read_value(src) {
        Ok(v) => v,
        Err(e) => return Some(Err(e)),
    };

    Some(Ok((key, Some(value))))
}

fn read_key<'a>(src: &mut &'a str) -> io::Result<(&'a str, bool)> {
    let s = src.as_bytes();

    let mut r#match = None;

    let key = if let Some(i) = memchr::memchr2(SEPARATOR, DELIMITER, s) {
        let (k, rest) = src.split_at(i);
        *src = &rest[1..];
        r#match = Some(s[i]);
        k
    } else {
        let (k, rest) = src.split_at(src.len());
        *src = rest;
        k
    };

    let is_delimited = matches!(r#match, Some(DELIMITER));

    if src.is_empty() && is_delimited {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unexpected field delimiter after key",
        ))
    } else if key.is_empty() {
        Err(io::Error::new(io::ErrorKind::InvalidData, "missing key"))
    } else {
        let is_separated = matches!(r#match, Some(SEPARATOR));
        Ok((key, is_separated))
    }
}

fn read_value<'a>(src: &mut &'a str) -> io::Result<&'a str> {
    let mut is_delimited = false;

    let value = if let Some(i) = memchr::memchr(DELIMITER, src.as_bytes()) {
        let (v, rest) = src.split_at(i);
        *src = &rest[1..];
        is_delimited = true;
        v
    } else {
        let (v, rest) = src.split_at(src.len());
        *src = rest;
        v
    };

    if src.is_empty() && is_delimited {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unexpected field delimiter after value",
        ))
    } else {
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() -> io::Result<()> {
        let mut src = "";
        assert!(next(&mut src).is_none());

        let mut src = "NS=2";
        assert_eq!(next(&mut src).transpose()?, Some(("NS", Some("2"))));

        let mut src = "NS=2;DP=.;H3";
        assert_eq!(next(&mut src).transpose()?, Some(("NS", Some("2"))));
        assert_eq!(next(&mut src).transpose()?, Some(("DP", Some("."))));
        assert_eq!(next(&mut src).transpose()?, Some(("H3", None)));

        // unexpected field delimiter after key
        let mut src = "H3;";
        assert!(matches!(
            next(&mut src),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData,
        ));

        // missing key
        let mut src = ";";
        assert!(matches!(
            next(&mut src),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData,
        ));

        // unexpected field delimiter after value
        let mut src = "NS=2;";
        assert!(matches!(
            next(&mut src),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
