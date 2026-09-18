use noodles_bam as bam;
use noodles_cram::crai;
use noodles_csi as csi;

/// An alignment index.
pub enum Index {
    /// A SAM index.
    Sam(csi::Index),
    /// A BAM index.
    Bam(bam::Index),
    /// A CRAM index.
    Cram(crai::Index),
}
