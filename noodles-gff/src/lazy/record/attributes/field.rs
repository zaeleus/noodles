//! Raw GFF record attributes field.

mod value;

use std::io;

use self::value::parse_value;
pub use self::value::Value;

pub(super) fn parse_field<'a>(buf: &mut &'a str) -> io::Result<(&'a str, Value<'a>)> {
    const DELIMITER: u8 = b';';
    const SEPARATOR: char = '=';

    let (raw_field, rest) = match buf.as_bytes().iter().position(|&b| b == DELIMITER) {
        Some(i) => {
            let (s, r) = buf.split_at(i);
            (s, &r[1..])
        }
        None => buf.split_at(buf.len()),
    };

    *buf = rest;

    let (key, raw_value) = raw_field
        .split_once(SEPARATOR)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid field"))?;

    let value = parse_value(raw_value);

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
