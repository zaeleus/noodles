//! Async FAI filesystem operations.

use std::path::Path;

use tokio::{
    fs::File,
    io::{self, BufReader, BufWriter},
};

use super::io::{Reader, Writer};
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

/// Writes a FASTA index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_fasta::fai;
/// let index = fai::Index::default();
/// fai::r#async::fs::write("reference.fa.fai", &index).await?;
/// # Ok(())
/// # }
/// ```
pub async fn write<P>(dst: P, index: &Index) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst)
        .await
        .map(BufWriter::new)
        .map(Writer::new)?;

    writer.write_index(index).await?;

    writer.shutdown().await?;

    Ok(())
}
