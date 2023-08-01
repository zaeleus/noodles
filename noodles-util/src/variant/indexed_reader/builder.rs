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
    ///
    /// By default, the compression is autodetected on build. This can be used to override it, but
    /// note that only bgzip-compressed streams can be indexed.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, Compression};
    /// let builder = variant::indexed_reader::Builder::default()
    ///     .set_compression(Some(Compression::Bgzf));
    /// ```
    pub fn set_compression(mut self, compression: Option<Compression>) -> Self {
        self.compression = Some(compression);
        self
    }

    /// Sets the format of the input.
    ///
    /// By default, the format is autodetected on build. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, Format};
    /// let builder = variant::indexed_reader::Builder::default().set_format(Format::Vcf);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets an index.
    ///
    /// When building from a path ([`Self::build_from_path`]), an associated index at `<src>.tbi`
    /// or `<src>.csi` will attempt to be loaded. This can be used to override it if the index
    /// cannot be found or when building from a reader (`[Self::build_from_reader]`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use noodles_util::variant;
    ///
    /// let index = csi::Index::default();
    /// let builder = variant::indexed_reader::Builder::default().set_index(index);
    /// ```
    pub fn set_index(mut self, index: csi::Index) -> Self {
        self.index = Some(index);
        self
    }

    /// Builds an indexed variant reader from a path.
    ///
    /// The compression method and format will be autodetected, if not overridden. If no index is
    /// set ([`Self::set_index`]), this will attempt to load an associated index at `<src>.tbi` or
    /// `<src>.csi`.
    ///
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_util::variant;
    /// let reader = variant::indexed_reader::Builder::default().build_from_path("sample.vcf.gz")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
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
    ///
    /// The compression method and format will be autodetected, if not overridden. An index must be
    /// set (`[Self::set_index]`). The reader must be a bgzip-compressed stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Write};
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi as csi;
    /// use noodles_util::variant;
    ///
    /// let mut writer = bgzf::Writer::new(Vec::new());
    /// writer.write_all(b"BCF")?;
    /// let data = writer.finish()?;
    ///
    /// let index = csi::Index::default();
    /// let reader = variant::indexed_reader::Builder::default()
    ///     .set_index(index)
    ///     .build_from_reader(&data[..])?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
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
