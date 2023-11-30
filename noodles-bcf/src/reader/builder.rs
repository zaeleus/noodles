use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Reader;

/// A BCF reader builder.
#[derive(Default)]
pub struct Builder;

impl Builder {
    /// Builds a BCF reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bcf::reader::Builder;
    /// let reader = Builder.build_from_path("sample.bcf")?;
    /// # Ok::<_, std::io::Error>(())
    /// ````
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<bgzf::Reader<File>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(src)?;
        self.build_from_reader(file)
    }

    /// Builds a BCF reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::reader::Builder;
    /// let reader = Builder::default().build_from_reader(io::empty());
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<Reader<bgzf::Reader<R>>>
    where
        R: Read,
    {
        Ok(Reader::new(reader))
    }
}
