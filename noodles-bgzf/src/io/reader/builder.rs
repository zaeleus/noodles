use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use super::Reader;
use crate::io::Block;

/// A BGZF reader builder.
#[derive(Debug, Default)]
pub struct Builder;

impl Builder {
    /// Builds a BGZF reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::io::reader::Builder::default().build_from_path("example.gz")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(src)?;
        Ok(self.build_from_reader(file))
    }

    /// Builds a BGZF reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::io::reader::Builder::default().build_from_reader(io::empty());
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: Read,
    {
        Reader {
            inner: reader,
            buf: Vec::new(),
            position: 0,
            block: Block::default(),
        }
    }
}
