//! CRAM index record and fields.

mod field;

pub use self::field::Field;

use std::{convert::TryFrom, error, fmt, str::FromStr};

use noodles_bam as bam;

const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 6;

/// A CRAM index record.
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
    /// Creates a CRAM index record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_bam as bam;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
    ///     10946,
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    /// # Ok::<(), bam::record::reference_sequence_id::TryFromIntError>(())
    /// ```
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

    /// Returns the reference sequence ID.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_bam as bam;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
    ///     10946,
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.reference_sequence_id().map(i32::from), Some(0));
    /// # Ok::<(), bam::record::reference_sequence_id::TryFromIntError>(())
    /// ```
    pub fn reference_sequence_id(&self) -> Option<bam::record::ReferenceSequenceId> {
        self.reference_sequence_id
    }

    /// Returns the alignment start.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_bam as bam;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
    ///     10946,
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.alignment_start(), 10946);
    /// # Ok::<(), bam::record::reference_sequence_id::TryFromIntError>(())
    /// ```
    pub fn alignment_start(&self) -> i32 {
        self.alignment_start
    }

    /// Returns the alignment span.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_bam as bam;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
    ///     10946,
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.alignment_span(), 6765);
    /// # Ok::<(), bam::record::reference_sequence_id::TryFromIntError>(())
    /// ```
    pub fn alignment_span(&self) -> i32 {
        self.alignment_span
    }

    /// Returns the offset of the container from the start of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_bam as bam;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
    ///     10946,
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.offset(), 17711);
    /// # Ok::<(), bam::record::reference_sequence_id::TryFromIntError>(())
    /// ```
    pub fn offset(&self) -> u64 {
        self.offset
    }

    /// Returns the offset of the slice from the start of the container.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_bam as bam;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
    ///     10946,
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.landmark(), 233);
    /// # Ok::<(), bam::record::reference_sequence_id::TryFromIntError>(())
    /// ```
    pub fn landmark(&self) -> u64 {
        self.landmark
    }

    /// Returns the size of the slice in bytes.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_bam as bam;
    /// use noodles_cram::crai;
    ///
    /// let record = crai::Record::new(
    ///     bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
    ///     10946,
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// assert_eq!(record.slice_length(), 317811);
    /// # Ok::<(), bam::record::reference_sequence_id::TryFromIntError>(())
    /// ```
    pub fn slice_length(&self) -> u64 {
        self.slice_length
    }
}

/// An error returned when a raw CRAM index record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    Missing(Field),
    /// A field is invalid.
    Invalid(Field, std::num::ParseIntError),
    /// The reference sequence ID is invalid.
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
        .ok_or(ParseError::Missing(field))
        .and_then(|s| s.parse().map_err(|e| ParseError::Invalid(field, e)))
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

    #[test]
    fn test_from_str_with_invalid_records() {
        assert_eq!(
            "0\t10946".parse::<Record>(),
            Err(ParseError::Missing(Field::AlignmentSpan))
        );

        assert!(matches!(
            "0\t10946\tnoodles".parse::<Record>(),
            Err(ParseError::Invalid(Field::AlignmentSpan, _))
        ));

        assert!(matches!(
            "-8\t10946\t6765\t17711\t233\t317811".parse::<Record>(),
            Err(ParseError::InvalidReferenceSequenceId(_))
        ));
    }
}
