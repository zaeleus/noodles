//! Async gzip index filesystem operations.

use std::path::Path;

use tokio::{
    fs::File,
    io::{self, AsyncWriteExt, BufReader, BufWriter},
};

use super::io::{Reader, Writer};
use crate::gzi::Index;

/// Reads the entire contents of a GZ index.
///
/// This is a convenience function and is equivalent to opening the given path and reading the
/// index.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_bgzf::gzi;
/// let index = gzi::r#async::fs::read("in.gz.gzi").await?;
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

/// Writes a GZ index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_bgzf::gzi;
/// let index = gzi::Index::default();
/// gzi::r#async::fs::write("in.gz.gzi", &index).await?;
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

    writer.get_mut().shutdown().await?;

    Ok(())
}
