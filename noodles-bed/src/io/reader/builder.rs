use std::{
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
};

use super::Reader;

/// A BED reader builder.
#[derive(Default)]
pub struct Builder<const N: usize>;

impl<const N: usize> Builder<N> {
    /// Builds a BED reader from a path.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::io::reader::Builder;
    /// let reader = Builder::<3>::default().build_from_path("in.bed");
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<N, BufReader<File>>>
    where
        P: AsRef<Path>,
    {
        File::open(src).map(BufReader::new).map(Reader::new)
    }

    /// Builds a BED reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed::io::reader::Builder;
    /// let reader = Builder::<3>::default().build_from_reader(io::empty());
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<N, BufReader<R>>
    where
        R: Read,
    {
        Reader::new(BufReader::new(reader))
    }
}
