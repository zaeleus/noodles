use std::io;

const SEPARATOR: char = ' ';
const DOUBLE_QUOTES: char = '"';
const TERMINATOR: char = ';';

pub(super) fn parse_field<'a>(s: &mut &'a str) -> io::Result<(&'a str, &'a str)> {
    let key = parse_key(s)?;
    let value = parse_value(s)?;
    maybe_consume_terminator(s);
    Ok((key, value))
}

fn parse_key<'a>(s: &mut &'a str) -> io::Result<&'a str> {
    let Some(i) = s.find(SEPARATOR) else {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    };

    let (t, rest) = s.split_at(i);
    *s = &rest[1..];
    Ok(t)
}

fn parse_value<'a>(s: &mut &'a str) -> io::Result<&'a str> {
    if let Some(rest) = s.strip_prefix(DOUBLE_QUOTES) {
        *s = rest;
        parse_string(s)
    } else {
        parse_raw_value(s)
    }
}

fn parse_string<'a>(s: &mut &'a str) -> io::Result<&'a str> {
    let Some(i) = s.find(DOUBLE_QUOTES) else {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    };

    let (t, rest) = s.split_at(i);
    *s = &rest[1..];
    Ok(t)
}

fn parse_raw_value<'a>(s: &mut &'a str) -> io::Result<&'a str> {
    let i = s.find(TERMINATOR).unwrap_or(s.len());
    let (t, rest) = s.split_at(i);
    *s = rest;
    Ok(t)
}

fn maybe_consume_terminator(s: &mut &str) {
    *s = s.trim_start();

    if let Some(rest) = s.strip_prefix(TERMINATOR) {
        *s = rest.trim_start();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() -> io::Result<()> {
        assert_eq!(parse_field(&mut r#"id "0""#)?, ("id", "0"));
        assert_eq!(parse_field(&mut "id 0")?, ("id", "0"));
        assert_eq!(parse_field(&mut r#"id "0";"#)?, ("id", "0"));
        assert_eq!(parse_field(&mut "id 0;")?, ("id", "0"));
        assert_eq!(parse_field(&mut r#"id "0" ; "#)?, ("id", "0"));
        assert_eq!(parse_field(&mut r#"id "0;1";"#)?, ("id", "0;1"));

        assert!(matches!(
            parse_field(&mut ""),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            parse_field(&mut "id"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
