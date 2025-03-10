use std::{error, fmt};

const SEPARATOR: char = ' ';
const DOUBLE_QUOTES: char = '"';
pub(super) const DELIMITER: char = ';';

/// An error returned when a raw GTF record attribute entry fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

pub(super) fn parse_field(s: &mut &str) -> Result<(String, String), ParseError> {
    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    let key = parse_key(s)?;
    let value = parse_value(s)?;
    discard_delimiter(s);

    Ok((key.into(), value.into()))
}

fn parse_key<'a>(s: &mut &'a str) -> Result<&'a str, ParseError> {
    let Some(i) = s.find(SEPARATOR) else {
        return Err(ParseError::Invalid);
    };

    let (key, rest) = s.split_at(i);
    *s = &rest[1..];

    Ok(key)
}

fn parse_value<'a>(s: &mut &'a str) -> Result<&'a str, ParseError> {
    if let Some(rest) = s.strip_prefix(DOUBLE_QUOTES) {
        *s = rest;
        parse_string(s)
    } else {
        parse_raw_value(s)
    }
}

fn parse_string<'a>(s: &mut &'a str) -> Result<&'a str, ParseError> {
    if let Some(i) = s.find(DOUBLE_QUOTES) {
        let (t, rest) = s.split_at(i);
        *s = &rest[1..];
        Ok(t)
    } else {
        Err(ParseError::Invalid)
    }
}

fn parse_raw_value<'a>(s: &mut &'a str) -> Result<&'a str, ParseError> {
    if let Some(i) = s.find(DELIMITER) {
        let (t, rest) = s.split_at(i);
        *s = rest;
        Ok(t)
    } else {
        Ok(s)
    }
}

fn discard_delimiter(s: &mut &str) {
    *s = s.trim_start();

    if let Some(rest) = s.strip_prefix(DELIMITER) {
        *s = rest;
    }

    *s = s.trim_start();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_field() {
        assert_eq!(
            parse_field(&mut r#"gene_id "g0""#),
            Ok((String::from("gene_id"), String::from("g0")))
        );
        assert_eq!(
            parse_field(&mut r#"gene_ids "g0;g1""#),
            Ok((String::from("gene_ids"), String::from("g0;g1")))
        );
        assert_eq!(
            parse_field(&mut r#"gene_id """#),
            Ok((String::from("gene_id"), String::from("")))
        );
        assert_eq!(
            parse_field(&mut r#"gene_id 0"#),
            Ok((String::from("gene_id"), String::from("0")))
        );

        assert_eq!(parse_field(&mut ""), Err(ParseError::Empty));
        assert_eq!(parse_field(&mut "gene_id"), Err(ParseError::Invalid));
        assert_eq!(parse_field(&mut r#""""#), Err(ParseError::Invalid));
    }
}
