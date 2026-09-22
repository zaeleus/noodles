//! BGZF filesystem operations.

use std::{io, path::Path};

use super::gzi;

/// Reads an associated gzip index.
///
/// This attempts to read an associated index at `<src>.gzi`.
///
/// # Examples
///
/// ```no_run
/// use noodles_bgzf as bgzf;
/// let index = bgzf::fs::read_associated_index("src.gz")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read_associated_index<P>(src: P) -> io::Result<gzi::Index>
where
    P: AsRef<Path>,
{
    const GZI_EXT: &str = "gzi";

    let gzi_src = src.as_ref().with_added_extension(GZI_EXT);
    gzi::fs::read(gzi_src)
}
