//! Raw GFF record attributes field.

mod tag;
mod value;

use std::{borrow::Cow, io, str};

use percent_encoding::percent_decode_str;

pub use self::value::Value;
use self::{tag::parse_tag, value::parse_value};

pub(super) fn parse_field<'a>(src: &mut &'a str) -> io::Result<(Cow<'a, str>, Value<'a>)> {
    const DELIMITER: char = ';';
    const SEPARATOR: char = '=';

    let (buf, rest) = src
        .split_once(DELIMITER)
        .unwrap_or_else(|| src.split_at(src.len()));

    *src = rest;

    if let Some((t, v)) = buf.split_once(SEPARATOR) {
        let tag = parse_tag(t)?;
        let value = parse_value(v)?;
        Ok((tag, value))
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidData, "invalid field"))
    }
}

fn percent_decode(s: &str) -> Result<Cow<'_, str>, str::Utf8Error> {
    percent_decode_str(s).decode_utf8()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() -> io::Result<()> {
        let mut src = "ID=1;Name=ndls";
        assert_eq!(
            parse_field(&mut src)?,
            (Cow::from("ID"), Value::String(Cow::from("1")))
        );
        assert_eq!(
            parse_field(&mut src)?,
            (Cow::from("Name"), Value::String(Cow::from("ndls")))
        );
        assert!(src.is_empty());

        let mut src = "ID";
        assert!(matches!(
            parse_field(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
