//! Async variant format filesystem operations.

use std::path::Path;

use noodles_bcf as bcf;
use noodles_vcf as vcf;
use tokio::io;

use crate::variant::{Index, fs::detect_format_from_extension, io::Format};

/// Reads an associated variant index.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_util::variant;
/// let index = variant::r#async::fs::read_associated_index("src.vcf.gz").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read_associated_index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    match detect_format_from_extension(&src) {
        Some(Format::Vcf) => vcf::r#async::fs::read_associated_index(src)
            .await
            .map(Index::Vcf),
        Some(Format::Bcf) => bcf::r#async::fs::read_associated_index(src)
            .await
            .map(Index::Bcf),
        None => Err(io::Error::from(io::ErrorKind::NotFound)),
    }
}
