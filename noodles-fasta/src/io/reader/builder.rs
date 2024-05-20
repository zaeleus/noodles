use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Reader;

/// A FASTA reader builder.
#[derive(Debug, Default)]
pub struct Builder;

impl Builder {
    /// Builds a FASTA reader from a path.
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<Box<dyn BufRead>>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        let file = File::open(src)?;

        let reader: Box<dyn BufRead> = match src.extension().and_then(|ext| ext.to_str()) {
            Some("gz" | "bgz") => Box::new(bgzf::Reader::new(file)),
            _ => Box::new(BufReader::new(file)),
        };

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
