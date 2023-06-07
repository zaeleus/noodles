use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use super::Reader;

/// A CRAM reader builder.
#[derive(Debug, Default)]
pub struct Builder;

impl Builder {
    /// Builds a CRAM reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_cram as cram;
    /// let reader = cram::reader::Builder::default().build_from_path("sample.cram")?;
    /// # Ok::<_, std::io::Error>(())
    /// ````
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        File::open(src).map(|file| self.build_from_reader(file))
    }

    /// Builds a CRAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    /// let reader = cram::reader::Builder::default().build_from_reader(io::empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: Read,
    {
        Reader::new(reader)
    }
}
