use noodles_csi as csi;
use noodles_vcf as vcf;

/// A variant index.
pub enum Index {
    /// A VCF index.
    Vcf(vcf::Index),
    /// A BCF index.
    Bcf(csi::Index),
}
