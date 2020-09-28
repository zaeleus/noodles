mod field;

pub use self::field::Field;

use std::{error, fmt, str::FromStr};

const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 6;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Record {
    reference_sequence_id: i32,
    alignment_start: i32,
    alignment_span: i32,
    offset: u64,
    landmark: u64,
    slice_len: u64,
}

impl Record {
    pub fn reference_sequence_id(&self) -> i32 {
        self.reference_sequence_id
    }

    pub fn alignment_start(&self) -> i32 {
        self.alignment_start
    }

    pub fn alignment_span(&self) -> i32 {
        self.alignment_span
    }

    pub fn offset(&self) -> u64 {
        self.offset
    }

    pub fn landmark(&self) -> u64 {
        self.landmark
    }

    pub fn slice_len(&self) -> u64 {
        self.slice_len
    }
}

#[derive(Debug)]
pub enum ParseError {
    Missing(Field),
    Invalid(Field, std::num::ParseIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing(field) => write!(f, "missing field: {:?}", field),
            Self::Invalid(field, message) => write!(f, "invalid {:?} field: {}", field, message),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

        let reference_sequence_id = parse_i32(&mut fields, Field::ReferenceSequenceId)?;
        let alignment_start = parse_i32(&mut fields, Field::AlignmentStart)?;
        let alignment_span = parse_i32(&mut fields, Field::AlignmentSpan)?;
        let offset = parse_u64(&mut fields, Field::Offset)?;
        let landmark = parse_u64(&mut fields, Field::Landmark)?;
        let slice_len = parse_u64(&mut fields, Field::SliceLength)?;

        Ok(Record {
            reference_sequence_id,
            alignment_start,
            alignment_span,
            offset,
            landmark,
            slice_len,
        })
    }
}

fn parse_i32<'a, I>(fields: &mut I, field: Field) -> Result<i32, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| ParseError::Missing(field))
        .and_then(|s| s.parse().map_err(|e| ParseError::Invalid(field, e)))
}

fn parse_u64<'a, I>(fields: &mut I, field: Field) -> Result<u64, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| ParseError::Missing(field))
        .and_then(|s| s.parse().map_err(|e| ParseError::Invalid(field, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Record = "0\t10946\t6765\t17711\t233\t317811".parse()?;

        let expected = Record {
            reference_sequence_id: 0,
            alignment_start: 10946,
            alignment_span: 6765,
            offset: 17711,
            landmark: 233,
            slice_len: 317811,
        };

        assert_eq!(actual, expected);

        Ok(())
    }
}
