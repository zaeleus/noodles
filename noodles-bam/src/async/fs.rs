//! Async BAM filesystem operations.

use std::path::Path;

use noodles_csi as csi;
use tokio::io;

use crate::{Index, bai, fs::path_with_added_extension};

/// Reads an associated BAM index.
///
/// This attempts to read an associated index at `<src>.bai` or `<src>.csi`, in that order.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_bam as bam;
/// let index = bam::r#async::fs::read_associated_index("src.bam").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read_associated_index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    const BAI_EXT: &str = "bai";
    const CSI_EXT: &str = "csi";

    let bai_src = path_with_added_extension(src.as_ref(), BAI_EXT);

    match bai::r#async::fs::read(bai_src).await {
        Ok(index) => Ok(Index::Bai(index)),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            let csi_src = path_with_added_extension(src.as_ref(), CSI_EXT);
            csi::r#async::fs::read(csi_src).await.map(Index::Csi)
        }
        Err(e) => Err(e),
    }
}
