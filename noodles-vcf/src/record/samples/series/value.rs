//! VCF record samples series value.

mod genotype;

use std::io;

pub use self::genotype::Genotype;
use crate::{
    io::reader::record_buf::value::percent_decode,
    variant::record::samples::series::{value::Array, Value},
    Header,
};

pub(crate) fn parse_value<'a>(
    src: &'a str,
    header: &Header,
    key: &str,
) -> io::Result<Option<Value<'a>>> {
    use crate::{
        header::record::value::map::format::{definition::definition, Number, Type},
        variant::record::samples::keys::key,
    };

    const MISSING: &str = ".";

    if src == MISSING {
        return Ok(None);
    } else if key == key::GENOTYPE {
        return parse_genotype_value(src).map(Some);
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
    let s = percent_decode(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if chars.next().is_none() {
            return Ok(Value::Character(c));
        }
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        "invalid character",
    ))
}

fn parse_string_value(src: &str) -> io::Result<Value<'_>> {
    percent_decode(src)
        .map(Value::String)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_genotype_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Genotype(Box::new(Genotype::new(src))))
}

fn parse_integer_array_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Array(Array::Integer(Box::new(src))))
}

fn parse_float_array_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Array(Array::Float(Box::new(src))))
}

fn parse_character_array_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Array(Array::Character(Box::new(src))))
}

fn parse_string_array_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::Array(Array::String(Box::new(src))))
}
