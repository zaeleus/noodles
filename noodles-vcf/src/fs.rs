//! VCF filesystem operations.

mod index;

pub use self::index::index;

use std::{
    ffi::OsStr,
    io,
    path::{Path, PathBuf},
};

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

    let tbi_src = path_with_added_extension(src.as_ref(), TABIX_EXT);

    match tabix::fs::read(tbi_src) {
        Ok(index) => Ok(Index::Tabix(index)),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            let csi_src = path_with_added_extension(src.as_ref(), CSI_EXT);
            csi::fs::read(csi_src).map(Index::Csi)
        }
        Err(e) => Err(e),
    }
}

pub(crate) fn path_with_added_extension<P, S>(src: P, ext: S) -> PathBuf
where
    P: AsRef<Path>,
    S: AsRef<OsStr>,
{
    let path = src.as_ref().to_path_buf();
    pathbuf_add_extension(path, ext)
}

fn pathbuf_add_extension<S>(src: PathBuf, ext: S) -> PathBuf
where
    S: AsRef<OsStr>,
{
    let mut s = src.into_os_string();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_with_added_extension() {
        assert_eq!(
            path_with_added_extension("sample.vcf.gz", "tbi"),
            PathBuf::from("sample.vcf.gz.tbi")
        );
    }
}
