//! Async alignment format filesystem operations.

use std::path::Path;

use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;
use tokio::io;

use crate::alignment::{Index, fs::detect_format_from_extension, io::Format};

/// Reads an associated alignment index.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_util::alignment;
/// let index = alignment::r#async::fs::read_associated_index("src.bam").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read_associated_index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    match detect_format_from_extension(&src) {
        Some(Format::Sam) => sam::r#async::fs::read_associated_index(src)
            .await
            .map(Index::Sam),
        Some(Format::Bam) => bam::r#async::fs::read_associated_index(src)
            .await
            .map(Index::Bam),
        Some(Format::Cram) => cram::r#async::fs::read_associated_index(src)
            .await
            .map(Index::Cram),
        None => Err(io::Error::from(io::ErrorKind::NotFound)),
    }
}
