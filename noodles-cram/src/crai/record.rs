mod field;

pub use self::field::Field;

use std::{convert::TryFrom, error, fmt, str::FromStr};

use noodles_bam as bam;

const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 6;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Record {
    reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
    alignment_start: i32,
    alignment_span: i32,
    offset: u64,
    landmark: u64,
    slice_length: u64,
}

impl Record {
    pub fn new(
        reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
        alignment_start: i32,
        alignment_span: i32,
        offset: u64,
        landmark: u64,
        slice_length: u64,
    ) -> Self {
        Self {
            reference_sequence_id,
            alignment_start,
            alignment_span,
            offset,
            landmark,
            slice_length,
        }
    }

    pub fn reference_sequence_id(&self) -> Option<bam::record::ReferenceSequenceId> {
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

    pub fn slice_length(&self) -> u64 {
        self.slice_length
    }
}

#[derive(Debug)]
pub enum ParseError {
    Missing(Field),
    Invalid(Field, std::num::ParseIntError),
    InvalidReferenceSequenceId(bam::record::reference_sequence_id::TryFromIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing(field) => write!(f, "missing field: {:?}", field),
            Self::Invalid(field, message) => write!(f, "invalid {:?} field: {}", field, message),
            Self::InvalidReferenceSequenceId(e) => {
                write!(f, "invalid reference sequence ID: {}", e)
            }
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

        let reference_sequence_id =
            parse_i32(&mut fields, Field::ReferenceSequenceId).and_then(|id| {
                if id == bam::record::reference_sequence_id::UNMAPPED {
                    Ok(None)
                } else {
                    bam::record::ReferenceSequenceId::try_from(id)
                        .map(Some)
                        .map_err(ParseError::InvalidReferenceSequenceId)
                }
            })?;

        let alignment_start = parse_i32(&mut fields, Field::AlignmentStart)?;
        let alignment_span = parse_i32(&mut fields, Field::AlignmentSpan)?;
        let offset = parse_u64(&mut fields, Field::Offset)?;
        let landmark = parse_u64(&mut fields, Field::Landmark)?;
        let slice_length = parse_u64(&mut fields, Field::SliceLength)?;

        Ok(Record {
            reference_sequence_id,
            alignment_start,
            alignment_span,
            offset,
            landmark,
            slice_length,
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
    fn test_from_str() -> Result<(), Box<dyn std::error::Error>> {
        let actual: Record = "0\t10946\t6765\t17711\t233\t317811".parse()?;

        let expected = Record {
            reference_sequence_id: Some(bam::record::ReferenceSequenceId::try_from(0)?),
            alignment_start: 10946,
            alignment_span: 6765,
            offset: 17711,
            landmark: 233,
            slice_length: 317811,
        };

        assert_eq!(actual, expected);

        Ok(())
    }
}
