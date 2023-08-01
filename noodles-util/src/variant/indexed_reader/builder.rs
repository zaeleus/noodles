use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_csi as csi;
use noodles_vcf as vcf;

use super::IndexedReader;
use crate::variant::{Compression, Format};

/// An indexed variant reader builder.
#[derive(Default)]
pub struct Builder {
    compression: Option<Option<Compression>>,
    format: Option<Format>,
    index: Option<csi::Index>,
}

impl Builder {
    /// Sets an index.
    pub fn set_index(mut self, index: csi::Index) -> Self {
        self.index = Some(index);
        self
    }

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
            (Format::Vcf, Some(Compression::Bgzf)) => {
                let mut builder = vcf::indexed_reader::Builder::default();

                if let Some(index) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_path(src).map(IndexedReader::Vcf)
            }
            (Format::Bcf, Some(Compression::Bgzf)) => {
                let mut builder = bcf::indexed_reader::Builder::default();

                if let Some(index) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_path(src).map(IndexedReader::Bcf)
            }
            (_, None) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "source not bgzip-compressed",
            )),
        }
    }
}
