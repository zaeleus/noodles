//! Variant format filesystem operations.

use std::{io, path::Path};

use noodles_bcf as bcf;
use noodles_vcf as vcf;

use super::{Index, io::Format};

/// Reads an associated variant index.
///
/// # Examples
///
/// ```no_run
/// use noodles_util::variant;
/// let index = variant::fs::read_associated_index("src.vcf.gz")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read_associated_index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    match detect_format_from_extension(&src) {
        Some(Format::Vcf) => vcf::fs::read_associated_index(src).map(Index::Vcf),
        Some(Format::Bcf) => bcf::fs::read_associated_index(src).map(Index::Bcf),
        None => Err(io::Error::from(io::ErrorKind::NotFound)),
    }
}

pub(crate) fn detect_format_from_extension<P>(src: P) -> Option<Format>
where
    P: AsRef<Path>,
{
    const VCF_GZ_SUFFIX: &str = ".vcf.gz";
    const BCF_SUFFIX: &str = ".bcf";

    src.as_ref()
        .file_name()
        .and_then(|filename| filename.to_str())
        .and_then(|filename| {
            if filename.ends_with(VCF_GZ_SUFFIX) {
                Some(Format::Vcf)
            } else if filename.ends_with(BCF_SUFFIX) {
                Some(Format::Bcf)
            } else {
                None
            }
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_format_from_extension() {
        assert_eq!(
            detect_format_from_extension("sample.vcf.gz"),
            Some(Format::Vcf)
        );
        assert_eq!(
            detect_format_from_extension("sample.bcf"),
            Some(Format::Bcf)
        );
        assert!(detect_format_from_extension("sample.txt").is_none());
    }
}
