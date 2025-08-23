//! Async BGZF filesystem operations.

use std::path::Path;

use tokio::{fs::File, io};

use super::io::Reader;

/// Opens a BGZF file.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_bgzf as bgzf;
/// let _reader = bgzf::r#async::fs::open("in.gz").await?;
/// # Ok(())
/// # }
/// ```
#[deprecated(
    since = "0.43.0",
    note = "Use `File::open(src).await.map(bgzf::async::io::Reader::new)` instead."
)]
pub async fn open<P>(src: P) -> io::Result<Reader<File>>
where
    P: AsRef<Path>,
{
    File::open(src).await.map(Reader::new)
}
