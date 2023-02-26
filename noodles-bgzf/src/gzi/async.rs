//! Async gzip index.

mod reader;

pub use self::reader::Reader;

use std::path::Path;

use tokio::{
    fs::File,
    io::{self, BufReader},
};

use super::Index;

/// Reads the entire contents of a GZ index.
///
/// This is a convenience function and is equivalent to opening the given path and reading the
/// index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// #
/// # #[tokio::main]
/// # async fn main() -> io::Result<()> {
/// use noodles_bgzf::gzi;
/// let index = gzi::r#async::read("in.gz.gzi").await?;
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
