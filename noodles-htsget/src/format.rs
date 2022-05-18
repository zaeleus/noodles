use serde::{Deserialize, Serialize};

/// A data format.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
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
