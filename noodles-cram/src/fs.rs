//! CRAM filesystem operations.

mod index;

use std::{io, path::Path};

pub use self::index::index;
use super::crai;

/// Reads an associated CRAM index.
///
/// This attempts to read an associated index at `<src>.crai`.
///
/// # Examples
///
/// ```no_run
/// use noodles_cram as cram;
/// let index = cram::fs::read_associated_index("src.cram")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read_associated_index<P>(src: P) -> io::Result<crai::Index>
where
    P: AsRef<Path>,
{
    const CRAI_EXT: &str = "crai";

    let crai_src = src.as_ref().with_added_extension(CRAI_EXT);
    crai::fs::read(crai_src)
}
