//! Async FAI filesystem operations.

use std::path::Path;

use tokio::{
    fs::File,
    io::{self, BufReader},
};

use super::io::Reader;
use crate::fai::Index;

/// Reads the entire contents of a FASTA index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// parsing each record.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_fasta::fai;
/// let index = fai::r#async::fs::read("reference.fa.fai").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).await.map(BufReader::new).map(Reader::new)?;
    reader.read_index().await
}
