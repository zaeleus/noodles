//! Alignment format filesystem operations.

use std::{io, path::Path};

use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;

use super::{Index, io::Format};

/// Reads an associated alignment index.
///
/// # Examples
///
/// ```no_run
/// use noodles_util::alignment;
/// let index = alignment::fs::read_associated_index("src.bam")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read_associated_index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    match detect_format_from_extension(&src) {
        Some(Format::Sam) => sam::fs::read_associated_index(src).map(Index::Sam),
        Some(Format::Bam) => bam::fs::read_associated_index(src).map(Index::Bam),
        Some(Format::Cram) => cram::fs::read_associated_index(src).map(Index::Cram),
        None => Err(io::Error::from(io::ErrorKind::NotFound)),
    }
}

pub(crate) fn detect_format_from_extension<P>(src: P) -> Option<Format>
where
    P: AsRef<Path>,
{
    const SAM_GZ_SUFFIX: &str = ".sam.gz";
    const BAM_SUFFIX: &str = ".bam";
    const CRAM_SUFFIX: &str = ".cram";

    src.as_ref()
        .file_name()
        .and_then(|filename| filename.to_str())
        .and_then(|filename| {
            if filename.ends_with(SAM_GZ_SUFFIX) {
                Some(Format::Sam)
            } else if filename.ends_with(BAM_SUFFIX) {
                Some(Format::Bam)
            } else if filename.ends_with(CRAM_SUFFIX) {
                Some(Format::Cram)
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
            detect_format_from_extension("sample.sam.gz"),
            Some(Format::Sam)
        );
        assert_eq!(
            detect_format_from_extension("sample.bam"),
            Some(Format::Bam)
        );
        assert_eq!(
            detect_format_from_extension("sample.cram"),
            Some(Format::Cram)
        );
        assert!(detect_format_from_extension("sample.txt").is_none());
    }
}
