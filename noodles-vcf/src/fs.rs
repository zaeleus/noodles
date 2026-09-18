//! VCF filesystem operations.

mod index;

pub use self::index::index;

use std::{io, path::Path};

use noodles_csi as csi;
use noodles_tabix as tabix;

use crate::Index;

/// Reads an associated VCF index.
///
/// This attempts to read an associated index at `<src>.tbi` or `<src>.csi`, in that order.
///
/// # Examples
///
/// ```no_run
/// use noodles_vcf as vcf;
/// let index = vcf::fs::read_associated_index("src.vcf.gz")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read_associated_index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    const TABIX_EXT: &str = "tbi";
    const CSI_EXT: &str = "csi";

    let src = src.as_ref();
    let tbi_src = src.with_added_extension(TABIX_EXT);

    match tabix::fs::read(tbi_src) {
        Ok(index) => Ok(Index::Tabix(index)),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            let csi_src = src.with_added_extension(CSI_EXT);
            csi::fs::read(csi_src).map(Index::Csi)
        }
        Err(e) => Err(e),
    }
}
