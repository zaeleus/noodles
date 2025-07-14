//! Variant format filesystem operations.

use std::{
    io::{self, BufRead},
    path::Path,
};

use super::io::Reader;

/// Opens a variant format file.
///
/// # Examples
///
/// ```no_run
/// use noodles_util::variant;
/// let reader = variant::fs::open("in.vcf")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn open<P>(src: P) -> io::Result<Reader<Box<dyn BufRead>>>
where
    P: AsRef<Path>,
{
    super::io::reader::Builder::default().build_from_path(src)
}
