mod kind;
mod parser;

use std::{error, fmt, str::FromStr};

use self::kind::Kind;

type Field = (String, String);

#[derive(Debug)]
pub enum Record {
    Info(Vec<Field>),
    Filter(Vec<Field>),
    Format(Vec<Field>),
}

#[derive(Debug)]
pub enum ParseError {
    Invalid,
    InvalidKind(kind::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid record: ")?;

        match self {
            Self::Invalid => f.write_str("could not parse record"),
            Self::InvalidKind(e) => write!(f, "invalid kind: {}", e),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (_, (key, values)) = parser::parse(s).map_err(|_| ParseError::Invalid)?;

        let kind = key.parse().map_err(ParseError::InvalidKind)?;

        match kind {
            Kind::Info => Ok(Self::Info(values)),
            Kind::Filter => Ok(Self::Filter(values)),
            Kind::Format => Ok(Self::Format(values)),
        }
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
