mod field;

use std::{error, fmt, io, str::FromStr};

use noodles_core::region::Interval;

use self::field::Field;

const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 5;

/// A FASTA index record.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Record {
    name: Vec<u8>,
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

impl Record {
    /// Creates a FASTA index record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let record = fai::Record::new("sq0", 8, 4, 80, 81);
    /// ```
    pub fn new<N>(name: N, length: u64, offset: u64, line_bases: u64, line_width: u64) -> Self
    where
        N: Into<Vec<u8>>,
    {
        Self {
            name: name.into(),
            length,
            offset,
            line_bases,
            line_width,
        }
    }

    /// Returns the record name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let record = fai::Record::new("sq0", 8, 4, 80, 81);
    /// assert_eq!(record.name(), b"sq0");
    /// ```
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Returns the length of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let record = fai::Record::new("sq0", 8, 4, 80, 81);
    /// assert_eq!(record.len(), 8);
    /// ```
    #[allow(clippy::len_without_is_empty)]
    #[deprecated(since = "0.23.0", note = "Use `Record::length` instead.")]
    pub fn len(&self) -> u64 {
        self.length()
    }

    /// Returns the length of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let record = fai::Record::new("sq0", 8, 4, 80, 81);
    /// assert_eq!(record.length(), 8);
    /// ```
    pub fn length(&self) -> u64 {
        self.length
    }

    /// Returns the offset from the start.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let record = fai::Record::new("sq0", 8, 4, 80, 81);
    /// assert_eq!(record.offset(), 4);
    /// ```
    pub fn offset(&self) -> u64 {
        self.offset
    }

    /// Returns the number of bases in a line.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let record = fai::Record::new("sq0", 8, 4, 80, 81);
    /// assert_eq!(record.line_bases(), 80);
    /// ```
    pub fn line_bases(&self) -> u64 {
        self.line_bases
    }

    /// Returns the number of characters in a line.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let record = fai::Record::new("sq0", 8, 4, 80, 81);
    /// assert_eq!(record.line_width(), 81);
    /// ```
    pub fn line_width(&self) -> u64 {
        self.line_width
    }

    /// Returns the start position of the given interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::{region::Interval, Position};
    /// use noodles_fasta::fai;
    ///
    /// let record = fai::Record::new("sq0", 10946, 4, 80, 81);
    /// let interval = Interval::from(..);
    ///
    /// assert_eq!(record.query(interval)?, 4);
    /// Ok::<_, std::io::Error>(())
    /// ```
    pub fn query(&self, interval: Interval) -> io::Result<u64> {
        let start = interval
            .start()
            .map(|position| usize::from(position) - 1)
            .unwrap_or_default();

        let start =
            u64::try_from(start).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let pos = self.offset()
            + start / self.line_bases() * self.line_width()
            + start % self.line_bases();

        Ok(pos)
    }
}

/// An error returned when a raw FASTA index record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A field is missing.
    MissingField(Field),
    /// A field is invalid.
    InvalidField(Field, std::num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::MissingField(field) => write!(f, "missing field: {field:?}"),
            Self::InvalidField(field, _) => write!(f, "invalid field: {field:?}"),
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
        let offset = parse_u64(&mut fields, Field::Offset)?;
        let line_bases = parse_u64(&mut fields, Field::LineBases)?;
        let line_width = parse_u64(&mut fields, Field::LineWidth)?;

        Ok(Self::new(name, len, offset, line_bases, line_width))
    }
}

fn parse_string<'a, I>(fields: &mut I, field: Field) -> Result<String, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingField(field))
        .map(|s| s.into())
}

fn parse_u64<'a, I>(fields: &mut I, field: Field) -> Result<u64, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingField(field))
        .and_then(|s| s.parse().map_err(|e| ParseError::InvalidField(field, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!(
            "sq0\t10946\t4\t80\t81".parse(),
            Ok(Record::new("sq0", 10946, 4, 80, 81))
        );

        assert_eq!("".parse::<Record>(), Err(ParseError::Empty));

        assert_eq!(
            "sq0".parse::<Record>(),
            Err(ParseError::MissingField(Field::Length))
        );

        assert!(matches!(
            "sq0\tnoodles".parse::<Record>(),
            Err(ParseError::InvalidField(Field::Length, _))
        ));
    }
}
