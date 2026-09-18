//! BAM filesystem operations.

mod index;

use std::{fs::File, io, path::Path};

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

    let src = src.as_ref();
    let bai_src = src.with_added_extension(BAI_EXT);

    match bai::fs::read(bai_src) {
        Ok(index) => Ok(Index::Bai(index)),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            let csi_src = src.with_added_extension(CSI_EXT);
            csi::fs::read(csi_src).map(Index::Csi)
        }
        Err(e) => Err(e),
    }
}

fn open<P>(src: P) -> io::Result<Reader<bgzf::io::Reader<File>>>
where
    P: AsRef<Path>,
{
    File::open(src).map(Reader::new)
}
