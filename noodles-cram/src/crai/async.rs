//! Async CRAM index.

pub mod fs;
pub mod io;

use std::path::Path;

use tokio::fs::File;

use self::io::Writer;
use super::{Index, Record};

/// Reads the entire contents of a CRAM index.
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
/// use noodles_cram::crai;
/// let index = crai::r#async::read("sample.cram.crai").await?;
/// # Ok(())
/// # }
/// ```
#[deprecated(since = "0.100.0", note = "Use `crai::r#async::fs::read` instead.")]
pub async fn read<P>(src: P) -> tokio::io::Result<Index>
where
    P: AsRef<Path>,
{
    fs::read(src).await
}

/// Writes a CRAM index to a file.
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
/// use noodles_cram::crai;
/// let index = crai::Index::default();
/// crai::r#async::write("sample.cram.crai", &index).await?;
/// # Ok(())
/// # }
/// ```
pub async fn write<P>(dst: P, index: &[Record]) -> tokio::io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).await.map(Writer::new)?;
    writer.write_index(index).await?;
    writer.shutdown().await?;
    Ok(())
}
