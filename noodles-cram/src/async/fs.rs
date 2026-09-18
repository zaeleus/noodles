//! Async CRAM filesystem operations.

use std::path::Path;

use tokio::io;

use crate::crai;

/// Reads an associated CRAM index.
///
/// This attempts to read an associated index at `<src>.crai`.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_cram as cram;
/// let index = cram::r#async::fs::read_associated_index("src.cram").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read_associated_index<P>(src: P) -> io::Result<crai::Index>
where
    P: AsRef<Path>,
{
    const CRAI_EXT: &str = "crai";

    let crai_src = src.as_ref().with_added_extension(CRAI_EXT);
    crai::r#async::fs::read(crai_src).await
}
