//! Async CSI.

pub mod io;

use std::path::Path;

use tokio::fs::File;

use super::Index;

#[deprecated(
    since = "0.39.0",
    note = "Use `noodles_csi::r#async::io::Reader` instead."
)]
pub use self::io::Reader;

#[deprecated(
    since = "0.39.0",
    note = "Use `noodles_csi::r#async::io::Writer` instead."
)]
pub use self::io::Writer;

/// Reads the entire contents of a coordinate-sorted index (CSI).
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// #
/// # #[tokio::main]
/// # async fn main() -> io::Result<()> {
/// use noodles_csi as csi;
/// let index = csi::r#async::read("sample.bcf.csi").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read<P>(src: P) -> tokio::io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).await.map(Reader::new)?;
    reader.read_index().await
}

/// Writes a coordinate-sorted index (CSI) to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// #
/// # #[tokio::main]
/// # async fn main() -> io::Result<()> {
/// use noodles_csi as csi;
/// let index = csi::Index::default();
/// csi::r#async::write("sample.bcf.csi", &index).await?;
/// # Ok(())
/// # }
/// ```
pub async fn write<P>(dst: P, index: &Index) -> tokio::io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).await.map(Writer::new)?;
    writer.write_index(index).await?;
    writer.shutdown().await?;
    Ok(())
}
