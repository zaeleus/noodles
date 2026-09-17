//! BAM filesystem operations.

mod index;

use std::{
    ffi::OsStr,
    fs::File,
    io,
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;
use noodles_csi as csi;

pub use self::index::index;
use super::io::Reader;
use crate::{Index, bai};

/// Reads an associated BAM index.
///
/// This attempts to read an associated index at `<src>.bai` or `<src>.csi`, in that order.
///
/// # Examples
///
/// ```no_run
/// use noodles_bam as bam;
/// let index = bam::fs::read_associated_index("src.bam")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read_associated_index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    const BAI_EXT: &str = "bai";
    const CSI_EXT: &str = "csi";

    let bai_src = path_with_added_extension(src.as_ref(), BAI_EXT);

    match bai::fs::read(bai_src) {
        Ok(index) => Ok(Index::Bai(index)),
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

fn open<P>(src: P) -> io::Result<Reader<bgzf::io::Reader<File>>>
where
    P: AsRef<Path>,
{
    File::open(src).map(Reader::new)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_with_added_extension() {
        assert_eq!(
            path_with_added_extension("sample.bam", "bai"),
            PathBuf::from("sample.bam.bai")
        );
    }
}
