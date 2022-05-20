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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_format_for_crate_format() {
        assert_eq!(crate::Format::from(Format::Bam), crate::Format::Bam);
        assert_eq!(crate::Format::from(Format::Cram), crate::Format::Cram);
    }
}
