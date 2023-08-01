use std::{
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_csi as csi;
use noodles_vcf as vcf;

use super::IndexedReader;
use crate::variant::{
    reader::builder::{detect_compression, detect_format},
    Compression, Format,
};

/// An indexed variant reader builder.
#[derive(Default)]
pub struct Builder {
    compression: Option<Option<Compression>>,
    format: Option<Format>,
    index: Option<csi::Index>,
}

impl Builder {
    /// Sets the compression of the input.
    pub fn set_compression(mut self, compression: Option<Compression>) -> Self {
        self.compression = Some(compression);
        self
    }

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

    /// Builds an indexed variant reader from a path.
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<IndexedReader<R>>
    where
        R: Read,
    {
        let mut buf_reader = BufReader::new(reader);

        let compression = match self.compression {
            Some(compression) => compression,
            None => detect_compression(&mut buf_reader)?,
        };

        let format = match self.format {
            Some(format) => format,
            None => detect_format(&mut buf_reader, compression)?,
        };

        let reader = buf_reader.into_inner();

        match (format, compression) {
            (Format::Vcf, Some(Compression::Bgzf)) => {
                let mut builder = vcf::indexed_reader::Builder::default();

                if let Some(index) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_reader(reader).map(IndexedReader::Vcf)
            }
            (Format::Bcf, Some(Compression::Bgzf)) => {
                let mut builder = bcf::indexed_reader::Builder::default();

                if let Some(index) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_reader(reader).map(IndexedReader::Bcf)
            }
            (_, None) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "source not bgzip-compressed",
            )),
        }
    }
}
