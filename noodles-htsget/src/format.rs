use serde::Deserialize;

/// A data format.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
#[serde(rename_all = "UPPERCASE")]
pub enum Format {
    /// BAM.
    Bam,
    /// CRAM.
    Cram,
    /// VCF.
    Vcf,
    /// BCF.
    Bcf,
}
