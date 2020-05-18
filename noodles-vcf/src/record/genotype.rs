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
