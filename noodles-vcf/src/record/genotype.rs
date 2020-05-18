pub mod field;

pub use self::field::Field;

use std::{error, fmt, ops::Deref};

use super::Format;

const DELIMITER: char = ':';

#[derive(Debug, Default, PartialEq)]
pub struct Genotype(Vec<Field>);

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid genotype: {}", self.0)
    }
}

impl Genotype {
    pub fn from_str_format(s: &str, format: &Format) -> Result<Self, ParseError> {
        s.split(DELIMITER)
            .zip(format.iter())
            .map(|(t, k)| Field::from_str_key(t, k))
            .collect::<Result<_, _>>()
            .map(Self)
            .map_err(|_| ParseError(s.into()))
    }
}

impl Deref for Genotype {
    type Target = [Field];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Genotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, field) in self.iter().enumerate() {
            if i > 0 {
                write!(f, "{}", DELIMITER)?
            }

            write!(f, "{}", field)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::record::format;

    use super::*;

    #[test]
    fn test_from_str_format() -> Result<(), Box<dyn std::error::Error>> {
        let format = "GT".parse()?;
        let actual = Genotype::from_str_format("0|0", &format)?;
        assert_eq!(actual.len(), 1);

        let format = "GT:GQ".parse()?;
        let actual = Genotype::from_str_format("0|0:13", &format)?;
        assert_eq!(actual.len(), 2);

        Ok(())
    }

    #[test]
    fn test_fmt() {
        let genotype = Genotype(vec![Field::new(
            format::Key::Genotype,
            field::Value::String(String::from("0|0")),
        )]);

        assert_eq!(genotype.to_string(), "0|0");

        let genotype = Genotype(vec![
            Field::new(
                format::Key::Genotype,
                field::Value::String(String::from("0|0")),
            ),
            Field::new(
                format::Key::ConditionalGenotypeQuality,
                field::Value::Integer(13),
            ),
        ]);

        assert_eq!(genotype.to_string(), "0|0:13");
    }
}
