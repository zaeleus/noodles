//! BCF filesystem operations.

mod index;

use std::{fs::File, io, path::Path};

use noodles_bgzf as bgzf;
use noodles_csi as csi;

pub use self::index::index;
use super::io::Reader;

/// Reads an associated BCF index.
///
/// This attempts to read an associated index at `<src>.csi`.
///
/// # Examples
///
/// ```no_run
/// use noodles_bcf as bcf;
/// let index = bcf::fs::read_associated_index("src.bcf")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read_associated_index<P>(src: P) -> io::Result<csi::Index>
where
    P: AsRef<Path>,
{
    const CSI_EXT: &str = "csi";

    let csi_src = src.as_ref().with_added_extension(CSI_EXT);
    csi::fs::read(csi_src)
}

fn open<P>(src: P) -> io::Result<Reader<bgzf::io::Reader<File>>>
where
    P: AsRef<Path>,
{
    File::open(src).map(Reader::new)
}
