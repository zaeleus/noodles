use std::{
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
};

use noodles_bam as bam;
use noodles_cram::{self as cram, crai};
use noodles_csi as csi;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::IndexedReader;
use crate::alignment::{
    reader::builder::{detect_compression_method, detect_format},
    CompressionMethod, Format,
};

/// An alignment index.
pub enum Index {
    /// CSI.
    Csi(csi::Index),
    /// CRAI.
    Crai(crai::Index),
}

impl From<csi::Index> for Index {
    fn from(index: csi::Index) -> Self {
        Self::Csi(index)
    }
}

impl From<crai::Index> for Index {
    fn from(index: crai::Index) -> Self {
        Self::Crai(index)
    }
}

/// An indexed alignment reader builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
    index: Option<Index>,
}

impl Builder {
    /// Sets the compression method.
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the input.
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets an index.
    pub fn set_index<I>(mut self, index: I) -> Self
    where
        I: Into<Index>,
    {
        self.index = Some(index.into());
        self
    }

    /// Builds an indexed alignment reader from a path.
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<File>>
    where
        P: AsRef<Path>,
    {
        let mut reader = File::open(src.as_ref()).map(BufReader::new)?;

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => detect_compression_method(&mut reader)?,
        };

        let format = match self.format {
            Some(format) => format,
            None => detect_format(&mut reader, compression_method)?,
        };

        match (format, compression_method) {
            (Format::Sam, None) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "source not bgzip-compressed",
            )),
            (Format::Sam, Some(CompressionMethod::Bgzf)) => sam::indexed_reader::Builder::default()
                .build_from_path(src)
                .map(IndexedReader::Sam),
            (Format::Bam, Some(CompressionMethod::Bgzf)) => bam::indexed_reader::Builder::default()
                .build_from_path(src)
                .map(IndexedReader::Bam),
            (Format::Cram, None) => cram::indexed_reader::Builder::default()
                .set_reference_sequence_repository(self.reference_sequence_repository)
                .build_from_path(src)
                .map(IndexedReader::Cram),
            (_, _) => todo!(),
        }
    }

    /// Builds an indexed alignment reader from a reader.
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<IndexedReader<BufReader<R>>>
    where
        R: Read,
    {
        let mut reader = BufReader::new(reader);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => detect_compression_method(&mut reader)?,
        };

        let format = match self.format {
            Some(format) => format,
            None => detect_format(&mut reader, compression_method)?,
        };

        match (format, compression_method) {
            (Format::Sam, None) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "source not bgzip-compressed",
            )),
            (Format::Sam, Some(CompressionMethod::Bgzf)) => {
                let mut builder = sam::indexed_reader::Builder::default();

                if let Some(Index::Csi(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_reader(reader).map(IndexedReader::Sam)
            }
            (Format::Bam, Some(CompressionMethod::Bgzf)) => {
                let mut builder = bam::indexed_reader::Builder::default();

                if let Some(Index::Csi(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_reader(reader).map(IndexedReader::Bam)
            }
            (Format::Cram, None) => {
                let mut builder = cram::indexed_reader::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository);

                if let Some(Index::Crai(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_reader(reader).map(IndexedReader::Cram)
            }
            (_, _) => todo!(),
        }
    }
}
