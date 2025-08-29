use std::path::Path;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, AsyncBufReadExt, AsyncRead, BufReader},
};

use super::Reader;
use crate::alignment::io::{CompressionMethod, Format};

/// An async alignment reader builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
}

impl Builder {
    /// Sets the compression method.
    ///
    /// By default, the compression method is autodetected on build. This can be used to override
    /// it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{r#async::io::reader::Builder, io::CompressionMethod};
    /// let _builder = Builder::default().set_compression_method(Some(CompressionMethod::Bgzf));
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
    /// use noodles_util::alignment::{r#async::io::reader::Builder, io::Format};
    /// let _builder = Builder::default().set_format(Format::Sam);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use noodles_util::alignment::{r#async::io::reader::Builder, io::Format};
    /// let repository = fasta::Repository::default();
    /// let _builder = Builder::default().set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Builds an async alignment reader from a path.
    ///
    /// By default, the format and compression method will be autodetected. This can be overridden
    /// by using [`Self::set_format`] and [`Self::set_compression_method`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::alignment::r#async::io::reader::Builder;
    /// let _reader = Builder::default().build_from_path("sample.bam").await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_path<P>(self, src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(src).await?;
        self.build_from_reader(file).await
    }

    /// Builds an async alignment reader from a reader.
    ///
    /// By default, the format and compression method will be autodetected. This can be overridden
    /// by using [`Self::set_format`] and [`Self::set_compression_method`].
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::alignment::r#async::io::reader::Builder;
    /// use tokio::io;
    /// let reader = Builder::default().build_from_reader(io::empty()).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_reader<R>(self, reader: R) -> io::Result<Reader<R>>
    where
        R: AsyncRead + Unpin + 'static,
    {
        use super::Inner;
        use crate::alignment::io::reader::builder::{detect_compression_method, detect_format};

        let mut reader = BufReader::new(reader);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => {
                let mut src = reader.fill_buf().await?;
                detect_compression_method(&mut src)?
            }
        };

        let format = match self.format {
            Some(format) => format,
            None => {
                let mut src = reader.fill_buf().await?;
                detect_format(&mut src, compression_method)?
            }
        };

        let inner = match (format, compression_method) {
            (Format::Sam, None) => Inner::Sam(sam::r#async::io::Reader::new(reader)),
            (Format::Sam, Some(CompressionMethod::Bgzf)) => Inner::SamGz(
                sam::r#async::io::Reader::new(bgzf::r#async::io::Reader::new(reader)),
            ),
            (Format::Bam, None) => Inner::BamRaw(bam::r#async::io::Reader::from(reader)),
            (Format::Bam, Some(CompressionMethod::Bgzf)) => {
                Inner::Bam(bam::r#async::io::Reader::new(reader))
            }
            (Format::Cram, None) => Inner::Cram(
                cram::r#async::io::reader::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository)
                    .build_from_reader(reader),
            ),
            (Format::Cram, Some(CompressionMethod::Bgzf)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid format compression method: CRAM cannot be bgzip-compressed",
                ));
            }
        };

        Ok(Reader(inner))
    }
}
