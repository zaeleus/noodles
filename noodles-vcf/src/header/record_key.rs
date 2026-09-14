use super::record::key::Other;

/// A key that identifies a VCF header record.
///
/// This is used to track the order in which records were added to a [`crate::Header`], which is
/// otherwise only retained within a category.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum RecordKey {
    /// An `ALT` record.
    AlternativeAllele(String),
    /// A `contig` record.
    Contig(String),
    /// A `FILTER` record.
    Filter(String),
    /// A `FORMAT` record.
    Format(String),
    /// An `INFO` record.
    Info(String),
    /// A nonstandard record.
    Other(Other),
}
