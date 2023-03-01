/// A variant format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    /// Variant Call Format (VCF).
    Vcf,
    /// BCF.
    Bcf,
}

/// A variant compression.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    /// BGZF compression.
    Bgzf,
}
