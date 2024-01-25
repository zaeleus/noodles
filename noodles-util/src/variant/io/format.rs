/// A variant format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    /// Variant Call Format (VCF).
    Vcf,
    /// BCF.
    Bcf,
}
