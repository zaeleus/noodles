pub mod field;

pub use self::field::Field;

use std::{error, fmt, str::FromStr};

const DELIMITER: char = '\t';

#[derive(Clone, Debug, Default)]
pub struct Data {
    fields: Vec<Field>,
}

impl Data {
    pub fn new(fields: Vec<Field>) -> Self {
        Self { fields }
    }

    pub fn fields(&self) -> &[Field] {
        &self.fields
    }

    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    pub fn len(&self) -> usize {
        self.fields.len()
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, field) in self.fields.iter().enumerate() {
            if i > 0 {
                f.write_str("\t")?;
            }

            write!(f, "{}", field)?;
        }

        Ok(())
    }
}

#[derive(Debug)]
pub struct ParseError(Box<dyn std::error::Error>);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl FromStr for Data {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(Self::default());
        }

        s.split(DELIMITER)
            .map(|t| t.parse().map_err(|e| ParseError(Box::new(e))))
            .collect::<Result<_, _>>()
            .map(Self::new)
    }
}

#[cfg(test)]
mod tests {
    use super::field::{Tag, Value};

    use super::*;

    #[test]
    fn test_fmt() {
        let data = Data::new(vec![
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
