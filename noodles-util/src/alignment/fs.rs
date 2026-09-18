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
    const SAM_EXT: &str = "sam";
    const BAM_EXT: &str = "bam";
    const CRAM_EXT: &str = "cram";

    src.as_ref().extension().and_then(|ext| match ext.to_str() {
        Some(SAM_EXT) => Some(Format::Sam),
        Some(BAM_EXT) => Some(Format::Bam),
        Some(CRAM_EXT) => Some(Format::Cram),
        _ => None,
    })
}
