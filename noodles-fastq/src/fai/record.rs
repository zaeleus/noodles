mod field;

use std::{error, fmt, str::FromStr};

use self::field::Field;

const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 6;

/// A FASTQ index record.
#[derive(Debug, Default, Eq, PartialEq)]
pub struct Record {
    name: String,
    length: u64,
    sequence_offset: u64,
    line_bases: u64,
    line_width: u64,
    quality_scores_offset: u64,
}

impl Record {
    /// Creates a FASTQ index record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let record = fai::Record::new("r0", 8, 4, 8, 9, 15);
    /// ```
    pub fn new<N>(
        name: N,
        length: u64,
        sequence_offset: u64,
        line_bases: u64,
        line_width: u64,
        quality_scores_offset: u64,
    ) -> Self
    where
        N: Into<String>,
    {
        Self {
            name: name.into(),
            length,
            sequence_offset,
            line_bases,
            line_width,
            quality_scores_offset,
        }
    }

    /// Returns the name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let record = fai::Record::new("r0", 8, 4, 8, 9, 15);
    /// assert_eq!(record.name(), "r0");
    /// ```
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the length of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let record = fai::Record::new("r0", 8, 4, 8, 9, 15);
    /// assert_eq!(record.len(), 8);
    /// ```
    #[allow(clippy::len_without_is_empty)]
    #[deprecated(since = "0.9.0", note = "Use `Record::length` instead.")]
    pub fn len(&self) -> u64 {
        self.length
    }

    /// Returns the length of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let record = fai::Record::new("r0", 8, 4, 8, 9, 15);
    /// assert_eq!(record.length(), 8);
    /// ```
    pub fn length(&self) -> u64 {
        self.length
    }

    /// Returns the offset to the sequence from the start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let record = fai::Record::new("r0", 8, 4, 8, 9, 15);
    /// assert_eq!(record.sequence_offset(), 4);
    /// ```
    pub fn sequence_offset(&self) -> u64 {
        self.sequence_offset
    }

    /// Returns the number of bases in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let record = fai::Record::new("r0", 8, 4, 8, 9, 15);
    /// assert_eq!(record.line_bases(), 8);
    /// ```
    pub fn line_bases(&self) -> u64 {
        self.line_bases
    }

    /// Returns the number of characters in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let record = fai::Record::new("r0", 8, 4, 8, 9, 15);
    /// assert_eq!(record.line_width(), 9);
    /// ```
    pub fn line_width(&self) -> u64 {
        self.line_width
    }

    /// Returns the offset to the quality scores from the start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let record = fai::Record::new("r0", 8, 4, 8, 9, 15);
    /// assert_eq!(record.quality_scores_offset(), 15);
    /// ```
    pub fn quality_scores_offset(&self) -> u64 {
        self.quality_scores_offset
    }
}

/// An error returned when a raw FASTQ index record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A field is missing.
    Missing(Field),
    /// A field is invalid.
    Invalid(Field, std::num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Missing(field) => write!(f, "missing field: {field:?}"),
            Self::Invalid(field, _) => write!(f, "invalid field: {field:?}"),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

        let name = parse_string(&mut fields, Field::Name)?;
        let len = parse_u64(&mut fields, Field::Length)?;
        let sequence_offset = parse_u64(&mut fields, Field::SequenceOffset)?;
        let line_bases = parse_u64(&mut fields, Field::LineBases)?;
        let line_width = parse_u64(&mut fields, Field::LineWidth)?;
        let quality_scores_offset = parse_u64(&mut fields, Field::QualityScoresOffset)?;

        Ok(Self {
            name,
            length: len,
            sequence_offset,
            line_bases,
            line_width,
            quality_scores_offset,
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
            "r0\t4\t4\t4\t5\t11".parse(),
            Ok(Record::new("r0", 4, 4, 4, 5, 11))
        );

        assert_eq!("".parse::<Record>(), Err(ParseError::Empty));

        assert_eq!(
            "r0".parse::<Record>(),
            Err(ParseError::Missing(Field::Length))
        );

        assert!(matches!(
            "r0\tnoodles".parse::<Record>(),
            Err(ParseError::Invalid(Field::Length, _))
        ));
    }
}
