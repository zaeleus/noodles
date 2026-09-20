//! Async BCF filesystem operations.

use std::path::Path;

use noodles_csi as csi;
use tokio::io;

/// Reads an associated BCF index.
///
/// This attempts to read an associated index at `<src>.csi`.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_bcf as bcf;
/// let index = bcf::r#async::fs::read_associated_index("src.bcf").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read_associated_index<P>(src: P) -> io::Result<csi::Index>
where
    P: AsRef<Path>,
{
    const CSI_EXT: &str = "csi";

    let csi_src = src.as_ref().with_added_extension(CSI_EXT);
    csi::r#async::fs::read(csi_src).await
}
