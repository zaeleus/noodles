//! VCF header record key.

use std::{error, fmt, str::FromStr};

/// A VCF header record key.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Key {
    /// File format (`fileformat`).
    FileFormat,
    /// Information (`INFO`).
    Info,
    /// Filter (`FILTER`)
    Filter,
    /// Genotype format (`FORMAT`).
    Format,
    /// Symbolic alternate allele (`ALT`).
    AlternativeAllele,
    /// Breakpoint assemblies URI (`assembly`).
    Assembly,
    /// Contig (`contig`).
    Contig,
    /// Any other record key.
    Other(String),
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Self::FileFormat => "fileformat",
            Self::Info => "INFO",
            Self::Filter => "FILTER",
            Self::Format => "FORMAT",
            Self::AlternativeAllele => "ALT",
            Self::Assembly => "assembly",
            Self::Contig => "contig",
            Self::Other(s) => s,
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw VCF header record kind fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
        }
    }
}

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "fileformat" => Ok(Self::FileFormat),
            "INFO" => Ok(Self::Info),
            "FILTER" => Ok(Self::Filter),
            "FORMAT" => Ok(Self::Format),
            "ALT" => Ok(Self::AlternativeAllele),
            "assembly" => Ok(Self::Assembly),
            "contig" => Ok(Self::Contig),
            _ => Ok(Self::Other(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Key::FileFormat.to_string(), "fileformat");
        assert_eq!(Key::Info.to_string(), "INFO");
        assert_eq!(Key::Filter.to_string(), "FILTER");
        assert_eq!(Key::Format.to_string(), "FORMAT");
        assert_eq!(Key::AlternativeAllele.to_string(), "ALT");
        assert_eq!(Key::Assembly.to_string(), "assembly");
        assert_eq!(Key::Contig.to_string(), "contig");
        assert_eq!(Key::Other(String::from("fileDate")).to_string(), "fileDate");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("fileformat".parse(), Ok(Key::FileFormat));
        assert_eq!("INFO".parse(), Ok(Key::Info));
        assert_eq!("FILTER".parse(), Ok(Key::Filter));
        assert_eq!("FORMAT".parse(), Ok(Key::Format));
        assert_eq!("ALT".parse(), Ok(Key::AlternativeAllele));
        assert_eq!("assembly".parse(), Ok(Key::Assembly));
        assert_eq!("contig".parse(), Ok(Key::Contig));
        assert_eq!("fileDate".parse(), Ok(Key::Other(String::from("fileDate"))));

        assert_eq!("".parse::<Key>(), Err(ParseError::Empty));
    }
}
