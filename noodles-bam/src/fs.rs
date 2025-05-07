//! BAM filesystem operations.

mod index;

use std::{fs::File, io, path::Path};

use noodles_bgzf as bgzf;

pub use self::index::index;
use super::io::Reader;

fn open<P>(src: P) -> io::Result<Reader<bgzf::io::Reader<File>>>
where
    P: AsRef<Path>,
{
    File::open(src).map(Reader::new)
}
