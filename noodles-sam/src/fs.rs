//! SAM filesystem operations.

mod index;

use std::{io, path::Path};

use noodles_csi as csi;

pub use self::index::index;

/// Reads an associated SAM index.
///
/// This attempts to read an associated index at `<src>.csi`.
///
/// # Examples
///
/// ```no_run
/// use noodles_sam as sam;
/// let index = sam::fs::read_associated_index("src.sam.gz")?;
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
