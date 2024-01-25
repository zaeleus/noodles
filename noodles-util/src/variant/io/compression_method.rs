/// A variant compression method.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CompressionMethod {
    /// BGZF compression.
    Bgzf,
}
