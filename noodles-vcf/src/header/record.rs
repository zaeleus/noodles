//! VCF header record and components.

pub mod key;
pub mod value;

pub use self::{key::Key, value::Value};

use std::str::FromStr;

use self::value::{
    map::{self, AlternativeAllele, Contig, Filter, Format, Info},
    Map,
};
use super::{parser::record::ParseError, FileFormat};

/// A VCF header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Record {
    /// An `ALT` record.
    AlternativeAllele(
        crate::record::alternate_bases::allele::Symbol,
        Map<AlternativeAllele>,
    ),
    /// A `contig` record.
    Contig(map::contig::Name, Map<Contig>),
    /// A `fileformat` record.
    FileFormat(FileFormat),
    /// A `FILTER` record.
    Filter(String, Map<Filter>),
    /// A `FORMAT` record.
    Format(crate::record::genotypes::keys::Key, Map<Format>),
    /// An `INFO` record.
    Info(crate::record::info::field::Key, Map<Info>),
    /// A nonstadard record.
    Other(key::Other, Value),
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_from((FileFormat::default(), s))
    }
}

impl TryFrom<(FileFormat, &str)> for Record {
    type Error = ParseError;

    fn try_from((file_format, s): (FileFormat, &str)) -> Result<Self, Self::Error> {
        super::parser::parse_record(s.as_bytes(), file_format)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let line = "##fileformat=VCFv4.3";
        assert_eq!(line.parse(), Ok(Record::FileFormat(FileFormat::new(4, 3))));

        let line =
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#;
        assert!(matches!(line.parse(), Ok(Record::Info(..))));

        assert!("".parse::<Record>().is_err());

        Ok(())
    }
}
