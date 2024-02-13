use std::io;

use super::Keys;
use crate::{
    variant::record::samples::series::{value::Array, Value},
    Header,
};

/// A VCF record samples sample.
#[derive(Debug, Eq, PartialEq)]
pub struct Sample<'a> {
    src: &'a str,
    keys: Keys<'a>,
}

impl<'a> Sample<'a> {
    pub(super) fn new(src: &'a str, keys: Keys<'a>) -> Self {
        Self { src, keys }
    }

    /// Returns the value at the given index.
    pub fn get_index<'h: 'a>(
        &self,
        header: &'h Header,
        i: usize,
    ) -> Option<Option<io::Result<Value<'a>>>> {
        self.values(header).nth(i)
    }

    pub fn values<'h: 'a>(
        &self,
        header: &'h Header,
    ) -> impl Iterator<Item = Option<io::Result<Value<'a>>>> + '_ {
        self.iter(header)
            .map(|result| result.map(|(_, value)| value).transpose())
    }

    /// Returns an iterator over fields.
    pub fn iter<'h: 'a>(
        &self,
        header: &'h Header,
    ) -> impl Iterator<Item = io::Result<(&str, Option<Value<'a>>)>> + '_ {
        const DELIMITER: char = ':';

        self.keys
            .iter()
            .zip(self.src.split(DELIMITER))
            .map(|(key, s)| parse_value(s, header, key).map(|value| (key, value)))
    }
}

impl<'a> AsRef<str> for Sample<'a> {
    fn as_ref(&self) -> &str {
        self.src
    }
}

fn parse_value<'a>(src: &'a str, header: &Header, key: &str) -> io::Result<Option<Value<'a>>> {
    use crate::header::{
        record::value::map::format::{definition::definition, Type},
        Number,
    };

    const MISSING: &str = ".";

    if src == MISSING {
        return Ok(None);
    }

    let (number, ty) = header
        .formats()
        .get(key)
        .map(|format| (format.number(), format.ty()))
        .or_else(|| definition(header.file_format(), key).map(|(n, t, _)| (n, t)))
        .unwrap_or_default();

    let value = match (number, ty) {
        (Number::Count(0), _) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid number for type",
            ))
        }
        (Number::Count(1), Type::Integer) => parse_integer_value(src)?,
        (Number::Count(1), Type::Float) => parse_float_value(src)?,
        (Number::Count(1), Type::Character) => parse_character_value(src)?,
        (Number::Count(1), Type::String) => parse_string_value(src)?,
        (_, Type::Integer) => parse_integer_array_value(src)?,
        (_, Type::Float) => parse_float_array_value(src)?,
        (_, Type::Character) => parse_character_array_value(src)?,
        (_, Type::String) => parse_string_array_value(src)?,
    };

    Ok(Some(value))
}

fn parse_integer_value(src: &str) -> io::Result<Value<'_>> {
    src.parse()
        .map(Value::Integer)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_float_value(src: &str) -> io::Result<Value<'_>> {
    src.parse()
        .map(Value::Float)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_character_value(src: &str) -> io::Result<Value<'_>> {
    let mut chars = src.chars();

    if let Some(c) = chars.next() {
        if chars.next().is_none() {
            return Ok(Value::Character(c));
        }
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        "invalid character value",
    ))
}

fn parse_string_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::String(src))
}

fn parse_integer_array_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Array(Array::Integer(Box::new(src))))
}

fn parse_float_array_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Array(Array::Integer(Box::new(src))))
}

fn parse_character_array_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Array(Array::Integer(Box::new(src))))
}

fn parse_string_array_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Array(Array::Integer(Box::new(src))))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_values() {
        let header = Header::default();
        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);
        let mut iter = sample.values(&header);

        assert!(matches!(
            iter.next(),
            Some(Some(Ok(Value::String(s)))) if s == "0|0"
        ));

        assert!(matches!(iter.next(), Some(None)));

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_iter() {
        let header = Header::default();
        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);
        let mut iter = sample.iter(&header);

        assert!(matches!(
            iter.next(),
            Some(Ok(("GT", Some(Value::String(s))))) if s == "0|0"
        ));

        assert!(matches!(iter.next(), Some(Ok(("GQ", None)))));

        assert!(iter.next().is_none());
    }
}
