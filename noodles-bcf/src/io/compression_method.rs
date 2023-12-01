/// A BCF compression method.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CompressionMethod {
    /// No compression.
    None,
    /// BGZF.
    Bgzf,
}
