/// A VCF index format.
#[derive(Clone, Copy, Default, Eq, PartialEq)]
pub enum Format {
    /// Coordinate-sorted index (CSI).
    Csi,
    /// tabix (TBI).
    #[default]
    Tabix,
}
