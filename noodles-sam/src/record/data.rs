//! SAM record data and fields.

pub mod field;

pub use self::field::Field;

use std::{error, fmt, str::FromStr};

const DELIMITER: char = '\t';

/// SAM record data.
///
/// This is also called optional fields.
#[derive(Clone, Debug, Default)]
pub struct Data {
    fields: Vec<Field>,
}

impl Data {
    /// Returns the list of data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let data = Data::from(vec![
    ///     Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
    ///     Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
    /// ]);
    ///
    /// let actual = data.fields();
    /// let expected = [
    ///     Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
    ///     Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
    /// ];
    /// assert_eq!(actual, expected);
    /// ```
    pub fn fields(&self) -> &[Field] {
        &self.fields
    }

    /// Returns whether there are any data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let data = Data::default();
    /// assert!(data.is_empty());
    ///
    /// let data = Data::from(vec![
    ///     Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
    ///     Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
    /// ]);
    /// assert!(!data.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    /// Returns the number of data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::{data::{field::{Tag, Value}, Field}, Data};
    ///
    /// let data = Data::default();
    /// assert_eq!(data.len(), 0);
    ///
    /// let data = Data::from(vec![
    ///     Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
    ///     Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
    /// ]);
    /// assert_eq!(data.len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        self.fields.len()
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, field) in self.fields.iter().enumerate() {
            if i > 0 {
                f.write_str("\t")?;
            }

            write!(f, "{}", field)?;
        }

        Ok(())
    }
}

impl From<Vec<Field>> for Data {
    fn from(fields: Vec<Field>) -> Self {
        Self { fields }
    }
}

/// An error returned when raw SAM record data fails to parse.
#[derive(Debug)]
pub enum ParseError {
    /// The input data contains an invalid field.
    InvalidField(field::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Data {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::default());
        }

        s.split(DELIMITER)
            .map(|t| t.parse().map_err(ParseError::InvalidField))
            .collect::<Result<Vec<_>, _>>()
            .map(Self::from)
    }
}

#[cfg(test)]
mod tests {
    use super::field::{Tag, Value};

    use super::*;

    #[test]
    fn test_fmt() {
        let data = Data::from(vec![
            Field::new(Tag::ReadGroup, Value::String(String::from("rg0"))),
            Field::new(Tag::AlignmentHitCount, Value::Int32(1)),
        ]);

        let expected = "RG:Z:rg0\tNH:i:1";

        assert_eq!(data.to_string(), expected);
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let fields = "RG:Z:rg0\tNH:i:1";
        let data: Data = fields.parse()?;
        let fields = data.fields();
        assert_eq!(fields.len(), 2);

        let fields = "";
        let data: Data = fields.parse()?;
        assert!(data.fields().is_empty());

        Ok(())
    }
}
