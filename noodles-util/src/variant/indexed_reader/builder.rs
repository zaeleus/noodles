use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_vcf as vcf;

use super::IndexedReader;
use crate::variant::{Compression, Format};

/// An indexed variant reader builder.
#[derive(Default)]
pub struct Builder {
    compression: Option<Option<Compression>>,
    format: Option<Format>,
}

impl Builder {
    /// Builds an indexed variant reader from a path.
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<File>>
    where
        P: AsRef<Path>,
    {
        use crate::variant::reader::builder::{detect_compression, detect_format};

        let mut reader = File::open(src.as_ref()).map(BufReader::new)?;

        let compression = match self.compression {
            Some(compression) => compression,
            None => detect_compression(&mut reader)?,
        };

        let format = match self.format {
            Some(format) => format,
            None => detect_format(&mut reader, compression)?,
        };

        match (format, compression) {
            (Format::Vcf, None) => todo!(),
            (Format::Vcf, Some(Compression::Bgzf)) => vcf::indexed_reader::Builder::default()
                .build_from_path(src)
                .map(IndexedReader::Vcf),
            (Format::Bcf, None) => todo!(),
            (Format::Bcf, Some(Compression::Bgzf)) => bcf::indexed_reader::Builder::default()
                .build_from_path(src)
                .map(IndexedReader::Bcf),
        }
    }
}
