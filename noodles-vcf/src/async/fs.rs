//! Async VCF filesystem operations.

use std::path::Path;

use noodles_csi as csi;
use noodles_tabix as tabix;

use tokio::io;

use crate::{Index, fs::path_with_added_extension};

/// Reads an associated VCF index.
///
/// This attempts to read an associated index at `<src>.tbi` or `<src>.csi`, in that order.
///
/// # Examples
///
/// ```no_run
/// # #[tokio::main]
/// # async fn main() -> tokio::io::Result<()> {
/// use noodles_vcf as vcf;
/// let index = vcf::r#async::fs::read_associated_index("src.vcf.gz").await?;
/// # Ok(())
/// # }
/// ```
pub async fn read_associated_index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    const TABIX_EXT: &str = "tbi";
    const CSI_EXT: &str = "csi";

    let tbi_src = path_with_added_extension(src.as_ref(), TABIX_EXT);

    match tabix::r#async::fs::read(tbi_src).await {
        Ok(index) => Ok(Index::Tabix(index)),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            let csi_src = path_with_added_extension(src.as_ref(), CSI_EXT);
            csi::r#async::fs::read(csi_src).await.map(Index::Csi)
        }
        Err(e) => Err(e),
    }
}
