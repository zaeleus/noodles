//! VCF record samples series value.

mod genotype;

use std::io;

pub use self::genotype::Genotype;
use crate::{
    io::reader::record_buf::value::percent_decode,
    variant::record::samples::{
        keys::key,
        series::{value::Array, Value},
    },
    Header,
};

pub(crate) fn parse_value<'a>(
    src: &'a str,
    header: &Header,
    key: &str,
) -> io::Result<Option<Value<'a>>> {
    use crate::header::record::value::map::format::{definition::definition, Number, Type};

    const MISSING: &str = ".";

    if src == MISSING {
        return Ok(None);
    } else if key == key::GENOTYPE {
        return Ok(Some(parse_genotype_value(src)));
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
        (_, Type::Integer) => parse_integer_array_value(src),
        (_, Type::Float) => parse_float_array_value(src),
        (_, Type::Character) => parse_character_array_value(src),
        (_, Type::String) => parse_string_array_value(src),
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

fn parse_genotype_value(src: &str) -> Value<'_> {
    Value::Genotype(Box::new(Genotype::new(src)))
}

fn parse_integer_array_value(src: &str) -> Value<'_> {
    Value::Array(Array::Integer(Box::new(src)))
}

fn parse_float_array_value(src: &str) -> Value<'_> {
    Value::Array(Array::Float(Box::new(src)))
}

fn parse_character_array_value(src: &str) -> Value<'_> {
    Value::Array(Array::Character(Box::new(src)))
}

fn parse_string_array_value(src: &str) -> Value<'_> {
    Value::Array(Array::String(Box::new(src)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        header::record::value::{
            map::{
                format::{Number, Type},
                Format,
            },
            Map,
        },
        variant::{
            record::samples::series::value::genotype::Phasing,
            record_buf::samples::sample::{value::genotype::Allele, Value as ValueBuf},
        },
    };

    #[test]
    fn test_parse_value() -> io::Result<()> {
        fn t(s: &str, header: &Header, key: &str, expected: Option<ValueBuf>) -> io::Result<()> {
            let actual = parse_value(s, header, key)
                .and_then(|result| result.map(ValueBuf::try_from).transpose())?;

            assert_eq!(actual, expected);

            Ok(())
        }

        #[rustfmt::skip]
        let header = Header::builder()
            .add_format("I32", Map::<Format>::new(Number::Count(1), Type::Integer, ""))
            .add_format("F32", Map::<Format>::new(Number::Count(1), Type::Float, ""))
            .add_format("CHAR", Map::<Format>::new(Number::Count(1), Type::Character, ""))
            .add_format("STRING", Map::<Format>::new(Number::Count(1), Type::String, ""))
            .add_format("I32_ARRAY", Map::<Format>::new(Number::Count(2), Type::Integer, ""))
            .add_format("F32_ARRAY", Map::<Format>::new(Number::Count(2), Type::Float, ""))
            .add_format("CHAR_ARRAY", Map::<Format>::new(Number::Count(2), Type::Character, ""))
            .add_format("STRING_ARRAY", Map::<Format>::new(Number::Count(2), Type::String, ""))
            .add_format("I32_INVALID", Map::<Format>::new(Number::Count(0), Type::Integer, ""))
            .build();

        t(".", &header, "I32", None)?;
        t("8", &header, "I32", Some(ValueBuf::from(8)))?;

        t(".", &header, "F32", None)?;
        t("0", &header, "F32", Some(ValueBuf::from(0.0)))?;

        t(".", &header, "CHAR", None)?;
        t("n", &header, "CHAR", Some(ValueBuf::from('n')))?;

        t(".", &header, "STRING", None)?;
        t("ndls", &header, "STRING", Some(ValueBuf::from("ndls")))?;

        t(".", &header, "I32_ARRAY", None)?;
        t(
            "8,.",
            &header,
            "I32_ARRAY",
            Some(ValueBuf::from(vec![Some(8), None])),
        )?;

        t(".", &header, "F32_ARRAY", None)?;
        t(
            "0,.",
            &header,
            "F32_ARRAY",
            Some(ValueBuf::from(vec![Some(0.0), None])),
        )?;

        t(".", &header, "CHAR_ARRAY", None)?;
        t(
            "n,.",
            &header,
            "CHAR_ARRAY",
            Some(ValueBuf::from(vec![Some('n'), None])),
        )?;

        t(".", &header, "STRING_ARRAY", None)?;
        t(
            "n,.",
            &header,
            "STRING_ARRAY",
            Some(ValueBuf::from(vec![Some(String::from("n")), None])),
        )?;

        t(".", &header, key::GENOTYPE, None)?;
        t(
            "0/0",
            &header,
            key::GENOTYPE,
            Some(ValueBuf::Genotype(
                [
                    Allele::new(Some(0), Phasing::Unphased),
                    Allele::new(Some(0), Phasing::Unphased),
                ]
                .into_iter()
                .collect(),
            )),
        )?;

        assert!(matches!(
            dbg!(parse_value("0", &header, "I32_INVALID")),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
