use std::{
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
};

use noodles_bam as bam;
use noodles_cram::{self as cram, crai};
use noodles_csi::{
    self as csi, BinningIndex,
    binning_index::index::reference_sequence::index::{BinnedIndex, LinearIndex},
};
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::IndexedReader;
use crate::alignment::io::{
    CompressionMethod, Format,
    reader::builder::{detect_compression_method, detect_format},
};

/// An alignment index.
pub enum Index {
    /// CSI.
    Csi(Box<dyn BinningIndex>),
    /// CRAI.
    Crai(crai::Index),
}

impl From<csi::binning_index::Index<BinnedIndex>> for Index {
    fn from(index: csi::binning_index::Index<BinnedIndex>) -> Self {
        Self::Csi(Box::new(index))
    }
}

impl From<csi::binning_index::Index<LinearIndex>> for Index {
    fn from(index: csi::binning_index::Index<LinearIndex>) -> Self {
        Self::Csi(Box::new(index))
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
    ///
    /// By default, the compression method is autodetected on build. This can be used to override
    /// it, but note that only bgzip-compressed streams can be indexed.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, io::CompressionMethod};
    /// let builder = alignment::io::indexed_reader::Builder::default()
    ///     .set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the input.
    ///
    /// By default, the format is autodetected on build. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, io::Format};
    /// let builder = alignment::io::indexed_reader::Builder::default()
    ///     .set_format(Format::Sam);
    /// ```
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
    ///
    /// When building from a path ([`Self::build_from_path`]), an associated index depending on the
    /// format will attempt to be loaded. This can be used to override it if the index cannot be
    /// found or when building from a reader ([`Self::build_from_reader`]).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use noodles_util::alignment;
    ///
    /// let index = bai::Index::default();
    /// let builder = alignment::io::indexed_reader::Builder::default()
    ///     .set_index(index);
    /// ```
    pub fn set_index<I>(mut self, index: I) -> Self
    where
        I: Into<Index>,
    {
        self.index = Some(index.into());
        self
    }

    /// Builds an indexed alignment reader from a path.
    ///
    /// The compression method and format will be autodetected, if not overridden. If no index is
    /// set ([`Self::set_index`]), this will attempt to load an associated index depending on the
    /// format.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_util::alignment;
    /// let reader = alignment::io::indexed_reader::Builder::default()
    ///     .build_from_path("sample.bam")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
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
            (Format::Sam | Format::Bam, None) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "source not bgzip-compressed",
            )),
            (Format::Sam, Some(CompressionMethod::Bgzf)) => {
                let mut builder = sam::io::indexed_reader::Builder::default();

                if let Some(Index::Csi(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_path(src).map(IndexedReader::Sam)
            }
            (Format::Bam, Some(CompressionMethod::Bgzf)) => {
                let mut builder = bam::io::indexed_reader::Builder::default();

                if let Some(Index::Csi(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_path(src).map(IndexedReader::Bam)
            }
            (Format::Cram, None) => {
                let mut builder = cram::io::indexed_reader::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository);

                if let Some(Index::Crai(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_path(src).map(IndexedReader::Cram)
            }
            (Format::Cram, Some(CompressionMethod::Bgzf)) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CRAM cannot be bgzip-compressed",
            )),
        }
    }

    /// Builds an indexed alignment reader from a reader.
    ///
    /// The compression method and format will be autodetected, if not overridden. An index must be
    /// set ([`Self::set_index`]). The reader must be a bgzip-compressed stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Write};
    /// use noodles_bam::bai;
    /// use noodles_bgzf as bgzf;
    /// use noodles_util::alignment;
    ///
    /// let mut writer = bgzf::io::Writer::new(Vec::new());
    /// writer.write_all(b"BAM\x01")?;
    /// let data = writer.finish()?;
    ///
    /// let index = bai::Index::default();
    /// let reader = alignment::io::indexed_reader::Builder::default()
    ///     .set_index(index)
    ///     .build_from_reader(&data[..])?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
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
            (Format::Sam | Format::Bam, None) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "source not bgzip-compressed",
            )),
            (Format::Sam, Some(CompressionMethod::Bgzf)) => {
                let mut builder = sam::io::indexed_reader::Builder::default();

                if let Some(Index::Csi(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_reader(reader).map(IndexedReader::Sam)
            }
            (Format::Bam, Some(CompressionMethod::Bgzf)) => {
                let mut builder = bam::io::indexed_reader::Builder::default();

                if let Some(Index::Csi(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_reader(reader).map(IndexedReader::Bam)
            }
            (Format::Cram, None) => {
                let mut builder = cram::io::indexed_reader::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository);

                if let Some(Index::Crai(index)) = self.index {
                    builder = builder.set_index(index);
                }

                builder.build_from_reader(reader).map(IndexedReader::Cram)
            }
            (Format::Cram, Some(CompressionMethod::Bgzf)) => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CRAM cannot be bgzip-compressed",
            )),
        }
    }
}
