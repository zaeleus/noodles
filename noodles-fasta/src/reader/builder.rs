use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use crate::Reader;

/// A FASTA reader builder.
#[derive(Debug, Default)]
pub struct Builder;

impl Builder {
    /// Builds a FASTA reader from a path.
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<BufReader<File>>>
    where
        P: AsRef<Path>,
    {
        let reader = File::open(&src).map(BufReader::new)?;
        self.build_from_reader(reader)
    }

    /// Builds a FASTA reader from a reader.
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<Reader<R>>
    where
        R: BufRead,
    {
        Ok(Reader::new(reader))
    }
}
