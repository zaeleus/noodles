/// The order in which VCF header records are written.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum RecordOrder {
    /// Records are grouped by category.
    ///
    /// Categories are written as `INFO`, `FILTER`, `FORMAT`, `ALT`, `contig`, and, finally,
    /// nonstandard records. Within a category, records keep the order they were added in.
    #[default]
    Grouped,
    /// Records are written in the order they were added to the header.
    ///
    /// This uses [`crate::Header::record_order`]. Records missing from it, e.g., those added using
    /// [`crate::Header::infos_mut`], are written last, grouped by category.
    AsAdded,
}
