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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_format_for_crate_format() {
        assert_eq!(crate::Format::from(Format::Vcf), crate::Format::Vcf);
        assert_eq!(crate::Format::from(Format::Bcf), crate::Format::Bcf);
    }
}
