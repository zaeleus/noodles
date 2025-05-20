//! Async tabix filesystem operations.

use std::path::Path;

use tokio::{fs::File, io};

use super::io::{Reader, Writer};
use crate::Index;

/// Reads the entire contents of a tabix index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> std::io::Result<()> {
/// use noodles_tabix as tabix;
/// let index = tabix::r#async::fs::read("sample.vcf.gz.tbi").await?;
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

/// Writes a tabix index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> std::io::Result<()> {
/// use noodles_csi::binning_index::index::Header;
/// use noodles_tabix as tabix;
///
/// let index = tabix::Index::builder().set_header(Header::default()).build();
/// tabix::r#async::fs::write("sample.vcf.gz.tbi", &index).await?;
/// # Ok(())
/// # }
/// ```
pub async fn write<P>(dst: P, index: &Index) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).await.map(Writer::new)?;
    writer.write_index(index).await?;
    writer.shutdown().await?;
    Ok(())
}
