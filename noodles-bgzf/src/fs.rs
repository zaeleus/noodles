//! BGZF filesystem operations.

use std::{fs::File, io, path::Path};

use super::io::Reader;

/// Opens a BGZF file.
///
/// # Examples
///
/// ```no_run
/// use noodles_bgzf as bgzf;
/// let _reader = bgzf::fs::open("in.gz")?;
/// # Ok::<_, std::io::Error>(())
/// ```
#[deprecated(
    since = "0.43.0",
    note = "Use `File::open(src).map(bgzf::io::Reader::new)` instead."
)]
pub fn open<P>(src: P) -> io::Result<Reader<File>>
where
    P: AsRef<Path>,
{
    File::open(src).map(Reader::new)
}
