/// A reads data format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    /// BAM.
    Bam,
    /// CRAM.
    Cram,
}

impl From<Format> for crate::Format {
    fn from(format: Format) -> Self {
        match format {
            Format::Bam => crate::Format::Bam,
            Format::Cram => crate::Format::Cram,
        }
    }
}
