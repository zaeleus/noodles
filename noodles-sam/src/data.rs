mod field;
mod subtype;
mod ty;
mod value;

pub use self::{field::Field, subtype::Subtype, ty::Type, value::Value};

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
    use super::*;

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
