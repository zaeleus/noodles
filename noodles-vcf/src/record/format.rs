use std::{error, fmt, ops::Deref, str::FromStr};

use crate::record::genotype::field::Key;

const DELIMITER: char = ':';

#[derive(Debug, Default)]
pub struct Format(Vec<Key>);

impl Deref for Format {
    type Target = [Key];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (i, key) in self.iter().enumerate() {
            if i > 0 {
                write!(f, "{}", DELIMITER)?
            }

            f.write_str(key.as_ref())?;
        }

        Ok(())
    }
}

#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid format: {}", self.0)
    }
}

impl FromStr for Format {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.split(DELIMITER)
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map(Format)
            .map_err(|_| ParseError(s.into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let format = Format(vec![Key::Genotype]);
        assert_eq!(format.to_string(), "GT");

        let format = Format(vec![
            Key::Genotype,
            Key::ConditionalGenotypeQuality,
            Key::ReadDepth,
            Key::HaplotypeQuality,
        ]);
        assert_eq!(format.to_string(), "GT:GQ:DP:HQ");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let actual: Format = "GT".parse()?;
        assert_eq!(actual.len(), 1);

        let actual: Format = "GT:GQ:DP:HQ".parse()?;
        assert_eq!(actual.len(), 4);

        Ok(())
    }
}
