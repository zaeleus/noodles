use std::io;

use crate::{
    header::{record::value::map::info::Type, Number},
    variant::record::info::field::{value::Array, Value},
};

pub(super) fn parse_value(src: &str, number: Number, ty: Type) -> io::Result<Value<'_>> {
    match (number, ty) {
        (Number::Count(0), Type::Flag) => parse_flag_value(src),
        (Number::Count(0), _) | (_, Type::Flag) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid number for type",
        )),
        (Number::Count(1), Type::Integer) => parse_integer_value(src),
        (Number::Count(1), Type::Float) => parse_float_value(src),
        (Number::Count(1), Type::Character) => parse_character_value(src),
        (Number::Count(1), Type::String) => parse_string_value(src),
        (_, Type::Integer) => parse_integer_array_value(src),
        (_, Type::Float) => parse_float_array_value(src),
        (_, Type::Character) => parse_character_array_value(src),
        (_, Type::String) => parse_string_array_value(src),
    }
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

fn parse_flag_value(src: &str) -> io::Result<Value<'_>> {
    if src.is_empty() {
        Ok(Value::Flag)
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidData, "invalid flag"))
    }
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
        "invalid character",
    ))
}

fn parse_string_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::String(src))
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
