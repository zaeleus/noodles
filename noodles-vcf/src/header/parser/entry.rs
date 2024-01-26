use crate::header::{
    record::value::{
        map::{AlternativeAllele, Contig, Filter, Format, Info},
        Map,
    },
    FileFormat,
};

/// A reference to an entry in the header.
pub enum Entry<'a> {
    /// A `fileformat` entry.
    FileFormat(FileFormat),
    /// An `INFO` entry.
    Info(&'a crate::record::info::field::Key, &'a Map<Info>),
    /// A `FILTER` entry.
    Filter(&'a str, &'a Map<Filter>),
    /// A `FORMAT` entry.
    Format(&'a crate::record::genotypes::keys::Key, &'a Map<Format>),
    /// An `ALT` entry.
    AlternativeAllele(
        &'a crate::record::alternate_bases::allele::Symbol,
        &'a Map<AlternativeAllele>,
    ),
    /// A `contig` entry.
    Contig(&'a str, &'a Map<Contig>),
    /// A nonstadard entry.
    Other,
    /// A header entry.
    Header,
}
