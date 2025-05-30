use crate::header::{
    FileFormat,
    record::value::{
        Map,
        map::{AlternativeAllele, Contig, Filter, Format, Info},
    },
};

/// A reference to an entry in the header.
pub enum Entry<'a> {
    /// A `fileformat` entry.
    FileFormat(FileFormat),
    /// An `INFO` entry.
    Info(&'a str, &'a Map<Info>),
    /// A `FILTER` entry.
    Filter(&'a str, &'a Map<Filter>),
    /// A `FORMAT` entry.
    Format(&'a str, &'a Map<Format>),
    /// An `ALT` entry.
    AlternativeAllele(&'a str, &'a Map<AlternativeAllele>),
    /// A `contig` entry.
    Contig(&'a str, &'a Map<Contig>),
    /// A nonstadard entry.
    Other,
    /// A header entry.
    Header,
}
