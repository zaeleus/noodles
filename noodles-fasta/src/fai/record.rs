mod field;

use std::{error, fmt, str::FromStr};

use self::field::Field;

const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 5;

/// A FASTA index record.
#[derive(Debug, Default, Eq, PartialEq)]
pub struct Record {
    reference_sequence_name: String,
    len: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

#[allow(clippy::len_without_is_empty)]
impl Record {
    /// Creates a FASTA index record.
    pub fn new(
        reference_sequence_name: String,
        len: u64,
        offset: u64,
        line_bases: u64,
        line_width: u64,
    ) -> Self {
        Self {
            reference_sequence_name,
            len,
            offset,
            line_bases,
            line_width,
        }
    }

    /// Returns the reference sequence name.
    pub fn reference_sequence_name(&self) -> &str {
        &self.reference_sequence_name
    }

    /// Returns the length of the sequence.
    pub fn len(&self) -> u64 {
        self.len
    }

    /// Returns the offset from the start.
    pub fn offset(&self) -> u64 {
        self.offset
    }

    /// Returns the number of bases in a line.
    pub fn line_bases(&self) -> u64 {
        self.line_bases
    }

    /// Returns the number of characters in a line.
    pub fn line_width(&self) -> u64 {
        self.line_width
    }
}

/// An error returned when a raw FASTA index record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    Missing(Field),
    /// A field is invalid.
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

        let reference_sequence_name = parse_string(&mut fields, Field::ReferenceSequenceName)?;
        let len = parse_u64(&mut fields, Field::Length)?;
        let offset = parse_u64(&mut fields, Field::Offset)?;
        let line_bases = parse_u64(&mut fields, Field::LineBases)?;
        let line_width = parse_u64(&mut fields, Field::LineWidth)?;

        Ok(Self {
            reference_sequence_name,
            len,
            offset,
            line_bases,
            line_width,
        })
    }
}

fn parse_string<'a, I>(fields: &mut I, field: Field) -> Result<String, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::Missing(field))
        .map(|s| s.into())
}

fn parse_u64<'a, I>(fields: &mut I, field: Field) -> Result<u64, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::Missing(field))
        .and_then(|s| s.parse().map_err(|e| ParseError::Invalid(field, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!(
            "sq0\t10946\t4\t80\t81".parse(),
            Ok(Record::new(String::from("sq0"), 10946, 4, 80, 81))
        );
    }
}
