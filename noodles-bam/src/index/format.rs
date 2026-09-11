/// A BAM index format.
#[derive(Clone, Copy, Default, Eq, PartialEq)]
pub enum Format {
    /// BAM index.
    #[default]
    Bai,
    /// Coordinate-sorted index (CSI).
    Csi,
}
