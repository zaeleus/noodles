//! VCF header record key.

pub mod other;

pub use self::other::Other;

use std::fmt;

/// VCF header record file format key.
pub const FILE_FORMAT: Key = Key::Standard(Standard::FileFormat);

/// VCF header record info key.
pub const INFO: Key = Key::Standard(Standard::Info);

/// VCF header record filter key.
pub const FILTER: Key = Key::Standard(Standard::Filter);

/// VCF header record format key.
pub const FORMAT: Key = Key::Standard(Standard::Format);

/// VCF header record alternative allele key.
pub const ALTERNATIVE_ALLELE: Key = Key::Standard(Standard::AlternativeAllele);

/// VCF header record assembly key.
pub const ASSEMBLY: Key = Key::Standard(Standard::Assembly);

/// VCF header record contig key.
pub const CONTIG: Key = Key::Standard(Standard::Contig);

/// VCF header record meta key.
pub const META: Key = Key::Standard(Standard::Meta);

/// A standard VCF record key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum Standard {
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
    /// Meta (`META`).
    Meta,
}

impl Standard {
    fn new(s: &str) -> Option<Self> {
        match s {
            "fileformat" => Some(Self::FileFormat),
            "INFO" => Some(Self::Info),
            "FILTER" => Some(Self::Filter),
            "FORMAT" => Some(Self::Format),
            "ALT" => Some(Self::AlternativeAllele),
            "assembly" => Some(Self::Assembly),
            "contig" => Some(Self::Contig),
            "META" => Some(Self::Meta),
            _ => None,
        }
    }
}

impl AsRef<str> for Standard {
    fn as_ref(&self) -> &str {
        match self {
            Self::FileFormat => "fileformat",
            Self::Info => "INFO",
            Self::Filter => "FILTER",
            Self::Format => "FORMAT",
            Self::AlternativeAllele => "ALT",
            Self::Assembly => "assembly",
            Self::Contig => "contig",
            Self::Meta => "META",
        }
    }
}

impl fmt::Display for Standard {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

/// A VCF header record key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum Key {
    /// A standard key.
    Standard(Standard),
    /// Any nonstandard key.
    Other(Other),
}

impl Key {
    /// Creates a nonstandard key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::header::record::Key;
    /// assert!(Key::other("INFO").is_none());
    /// assert!(Key::other("fileDate").is_some());
    /// ```
    pub fn other(s: &str) -> Option<Other> {
        match Self::from(s) {
            Self::Standard(_) => None,
            Self::Other(tag) => Some(tag),
        }
    }
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Self::Standard(tag) => tag.as_ref(),
            Self::Other(tag) => tag.as_ref(),
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

impl From<&str> for Key {
    fn from(s: &str) -> Self {
        match Standard::new(s) {
            Some(tag) => Self::Standard(tag),
            None => Self::Other(Other(s.into())),
        }
    }
}

impl From<String> for Key {
    fn from(s: String) -> Self {
        match Standard::new(&s) {
            Some(tag) => Self::Standard(tag),
            None => Self::Other(Other(s)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(FILE_FORMAT.to_string(), "fileformat");
        assert_eq!(INFO.to_string(), "INFO");
        assert_eq!(FILTER.to_string(), "FILTER");
        assert_eq!(FORMAT.to_string(), "FORMAT");
        assert_eq!(ALTERNATIVE_ALLELE.to_string(), "ALT");
        assert_eq!(ASSEMBLY.to_string(), "assembly");
        assert_eq!(CONTIG.to_string(), "contig");
        assert_eq!(META.to_string(), "META");
        assert_eq!(
            Key::Other(Other(String::from("fileDate"))).to_string(),
            "fileDate"
        );
    }

    #[test]
    fn test_from() {
        assert_eq!(Key::from("fileformat"), FILE_FORMAT);
        assert_eq!(Key::from("INFO"), INFO);
        assert_eq!(Key::from("FILTER"), FILTER);
        assert_eq!(Key::from("FORMAT"), FORMAT);
        assert_eq!(Key::from("ALT"), ALTERNATIVE_ALLELE);
        assert_eq!(Key::from("assembly"), ASSEMBLY);
        assert_eq!(Key::from("contig"), CONTIG);
        assert_eq!(Key::from("META"), META);
        assert_eq!(
            Key::from("fileDate"),
            Key::Other(Other(String::from("fileDate")))
        );
    }
}
