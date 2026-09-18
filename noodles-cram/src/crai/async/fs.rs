//! Async CRAM index filesystem operations.

use std::path::Path;

use tokio::{fs::File, io};

use super::io::Reader;
use crate::crai::Index;

/// Reads the entire contents of a CRAM index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_cram::crai;
/// let index = crai::r#async::fs::read("sample.cram.crai").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).await.map(Reader::new)?;
    reader.read_index().await
}
