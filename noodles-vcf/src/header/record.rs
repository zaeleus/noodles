mod kind;
mod parser;
mod value;

pub use self::kind::Kind;

use std::{error, fmt, str::FromStr};

use self::value::Value;

type Field = (String, String);

#[derive(Debug)]
pub enum Record {
    FileFormat(String),
    Info(Vec<Field>),
    Filter(Vec<Field>),
    Format(Vec<Field>),
    AlternativeAllele(Vec<Field>),
    Assembly(String),
    Contig(Vec<Field>),
    Other(String, String),
}

#[derive(Debug)]
pub enum ParseError {
    Invalid,
    InvalidKind(kind::ParseError),
    InvalidType,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid record: ")?;

        match self {
            Self::Invalid => f.write_str("could not parse record"),
            Self::InvalidKind(e) => write!(f, "invalid kind: {}", e),
            Self::InvalidType => f.write_str("invalid type"),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (_, (key, value)) = parser::parse(s).map_err(|_| ParseError::Invalid)?;

        let kind = key.parse().map_err(ParseError::InvalidKind)?;

        match kind {
            Kind::FileFormat => Ok(parse_string(value).map(Self::FileFormat)?),
            Kind::Info => Ok(parse_struct(value).map(Self::Info)?),
            Kind::Filter => Ok(parse_struct(value).map(Self::Filter)?),
            Kind::Format => Ok(parse_struct(value).map(Self::Format)?),
            Kind::AlternativeAllele => Ok(parse_struct(value).map(Self::AlternativeAllele)?),
            Kind::Assembly => Ok(parse_string(value).map(Self::Assembly)?),
            Kind::Contig => Ok(parse_struct(value).map(Self::Contig)?),
            Kind::Other(k) => {
                let v = parse_string(value)?;
                Ok(Self::Other(k, v))
            }
        }
    }
}

fn parse_string(value: Value) -> Result<String, ParseError> {
    match value {
        Value::String(v) => Ok(v),
        _ => Err(ParseError::InvalidType),
    }
}

fn parse_struct(value: Value) -> Result<Vec<(String, String)>, ParseError> {
    match value {
        Value::Struct(v) => Ok(v),
        _ => Err(ParseError::InvalidType),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let line =
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">"#;
        let record = line.parse()?;
        let expected = vec![
            (String::from("ID"), String::from("NS")),
            (String::from("Number"), String::from("1")),
            (String::from("Type"), String::from("Integer")),
            (
                String::from("Description"),
                String::from("Number of Samples With Data"),
            ),
        ];
        assert!(matches!(record, Record::Info(actual) if actual == expected));

        assert!("".parse::<Record>().is_err());

        Ok(())
    }
}
