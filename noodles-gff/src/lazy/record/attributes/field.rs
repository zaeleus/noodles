//! Raw GFF record attributes field.

mod value;

use std::io;

use self::value::parse_value;
pub use self::value::Value;

pub(super) fn parse_field<'a>(src: &mut &'a str) -> io::Result<(&'a str, Value<'a>)> {
    const DELIMITER: char = ';';
    const SEPARATOR: char = '=';

    let (buf, rest) = src
        .split_once(DELIMITER)
        .unwrap_or_else(|| src.split_at(src.len()));

    *src = rest;

    let (key, value) = buf
        .split_once(SEPARATOR)
        .map(|(k, v)| (k, parse_value(v)))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid field"))?;

    Ok((key, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() -> io::Result<()> {
        let mut src = "ID=nd;Name=ls";
        assert_eq!(parse_field(&mut src)?, ("ID", Value::String("nd")));
        assert_eq!(parse_field(&mut src)?, ("Name", Value::String("ls")));
        assert!(src.is_empty());

        let mut src = "ID";
        assert!(matches!(
            parse_field(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
