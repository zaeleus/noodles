mod kind;

use std::{error, fmt, str::FromStr};

pub use self::kind::Kind;

type Field = (String, String);

const DELIMITER: char = '\t';
const DATA_FIELD_DELIMITER: char = ':';

#[derive(Debug, Eq, PartialEq)]
pub enum Record {
    Header(Vec<Field>),
    ReferenceSequence(Vec<Field>),
    ReadGroup(Vec<Field>),
    Program(Vec<Field>),
    Comment(String),
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid record: {}", self.0)
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut pieces = s.split(DELIMITER);

        let kind = pieces
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| ParseError(s.into()))?;

        let record = if let Kind::Comment = kind {
            let comment = pieces.next().unwrap();
            Self::Comment(comment.into())
        } else {
            let fields = pieces
                .map(|field| {
                    let mut field_pieces = field.splitn(2, DATA_FIELD_DELIMITER);
                    let tag = field_pieces.next().ok_or_else(|| ParseError(s.into()))?;
                    let value = field_pieces.next().ok_or_else(|| ParseError(s.into()))?;
                    Ok((tag.into(), value.into()))
                })
                .collect::<Result<_, _>>()?;

            match kind {
                Kind::Header => Self::Header(fields),
                Kind::ReferenceSequence => Self::ReferenceSequence(fields),
                Kind::ReadGroup => Self::ReadGroup(fields),
                Kind::Program => Self::Program(fields),
                Kind::Comment => unreachable!(),
            }
        };

        Ok(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let record = "@HD\tVN:1.6".parse()?;
        let fields = vec![(String::from("VN"), String::from("1.6"))];
        assert!(matches!(record, Record::Header(f) if f == fields));

        let record = "@SQ\tSN:sq0\tLN:8".parse()?;
        let fields = vec![
            (String::from("SN"), String::from("sq0")),
            (String::from("LN"), String::from("8")),
        ];
        assert!(matches!(record, Record::ReferenceSequence(f) if f == fields));

        let record = "@RG\tID:rg0".parse()?;
        let fields = vec![(String::from("ID"), String::from("rg0"))];
        assert!(matches!(record, Record::ReadGroup(f) if f == fields));

        let record = "@PG\tID:pg0".parse()?;
        let fields = vec![(String::from("ID"), String::from("pg0"))];
        assert!(matches!(record, Record::Program(f) if f == fields));

        let record = "@CO\tnoodles".parse()?;
        assert!(matches!(record, Record::Comment(c) if c == "noodles"));

        Ok(())
    }
}
