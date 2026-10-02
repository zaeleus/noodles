//! Async FASTA filesystem operations.

use std::{io, path::Path};

use crate::fai;

/// Reads an associated FASTA index.
///
/// This attempts to read an associated index at `<src>.fai`.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_fasta as fasta;
/// let index = fasta::r#async::fs::read_associated_index("src.fa").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read_associated_index<P>(src: P) -> io::Result<fai::Index>
where
    P: AsRef<Path>,
{
    const FAI_EXT: &str = "fai";

    let fai_src = src.as_ref().with_added_extension(FAI_EXT);
    fai::r#async::fs::read(fai_src).await
}
