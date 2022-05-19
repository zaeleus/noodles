/// A variants data format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    /// VCF.
    Vcf,
    /// BCF.
    Bcf,
}

impl From<Format> for crate::Format {
    fn from(format: Format) -> Self {
        match format {
            Format::Vcf => crate::Format::Vcf,
            Format::Bcf => crate::Format::Bcf,
        }
    }
}
