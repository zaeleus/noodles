use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Reader;

/// A BAM reader builder.
#[derive(Debug, Default)]
pub struct Builder;

impl Builder {
    /// Builds a BAM reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_bam::io::reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.bam")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<bgzf::io::Reader<File>>>
    where
        P: AsRef<Path>,
    {
        File::open(src).map(|file| self.build_from_reader(file))
    }

    /// Builds a BAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::io::reader::Builder;
    /// let reader = Builder::default().build_from_reader(io::empty());
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<bgzf::io::Reader<R>>
    where
        R: Read,
    {
        Reader::new(reader)
    }
}
