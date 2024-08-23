use std::path::Path;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, AsyncWrite},
};

use super::Writer;
use crate::alignment::io::{CompressionMethod, Format};

/// An async alignment writer builder.
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
    /// use noodles_util::alignment::{r#async::io::writer::Builder, io::CompressionMethod};
    /// let builder = Builder::default().set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the output.
    ///
    /// By default, the format is autodetected on build. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{r#async::io::writer::Builder, io::Format};
    /// let builder = Builder::default().set_format(Format::Sam);
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
    /// use noodles_util::alignment::{r#async::io::writer::Builder, io::Format};
    /// let repository = fasta::Repository::default();
    /// let builder = Builder::default().set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Builds an async alignment writer from a path.
    ///
    /// By default, the format and compression method will be detected from the path extension.
    /// This can be overridden by using [`Self::set_format`] and [`Self::set_compression_method`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::alignment::r#async::io::reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.bam").await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_path<P>(
        mut self,
        src: P,
    ) -> io::Result<Writer<Box<dyn AsyncWrite + Unpin>>>
    where
        P: AsRef<Path>,
    {
        use crate::alignment::io::writer::builder::{
            detect_compression_method_from_path_extension, detect_format_from_path_extension,
        };

        let src = src.as_ref();

        if self.compression_method.is_none() {
            self.compression_method = Some(detect_compression_method_from_path_extension(src));
        }

        if self.format.is_none() {
            self.format = detect_format_from_path_extension(src);
        }

        File::create(src)
            .await
            .map(|file| self.build_from_writer(file))?
            .await
    }

    /// Builds an async alignment writer from a writer.
    ///
    /// If the format is not set, a default format is used. If the compression method is not set, a
    /// default one is determined by the format.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::alignment::r#async::io::writer::Builder;
    /// use tokio::io;
    ///
    /// let reader = Builder::default().build_from_writer(io::sink()).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_writer<W>(
        self,
        writer: W,
    ) -> io::Result<Writer<Box<dyn AsyncWrite + Unpin>>>
    where
        W: AsyncWrite + Unpin + 'static,
    {
        use bam::r#async::io::Writer as AsyncBamWriter;
        use sam::r#async::io::Writer as AsyncSamWriter;

        let format = self.format.unwrap_or(Format::Sam);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => match format {
                Format::Sam | Format::Cram => None,
                Format::Bam => Some(CompressionMethod::Bgzf),
            },
        };

        let writer = match (format, compression_method) {
            (Format::Sam, None) => {
                let inner: Box<dyn AsyncWrite + Unpin> = Box::new(writer);
                Writer::Sam(AsyncSamWriter::new(inner))
            }
            (Format::Sam, Some(CompressionMethod::Bgzf)) => {
                let encoder: Box<dyn AsyncWrite + Unpin> = Box::new(bgzf::AsyncWriter::new(writer));
                Writer::Sam(AsyncSamWriter::new(encoder))
            }
            (Format::Bam, None) => {
                let inner: Box<dyn AsyncWrite + Unpin> = Box::new(writer);
                Writer::Bam(AsyncBamWriter::from(inner))
            }
            (Format::Bam, Some(CompressionMethod::Bgzf)) => {
                let encoder: Box<dyn AsyncWrite + Unpin> = Box::new(bgzf::AsyncWriter::new(writer));
                Writer::Bam(AsyncBamWriter::from(encoder))
            }
            (Format::Cram, None) => {
                let inner: Box<dyn AsyncWrite + Unpin> = Box::new(writer);
                let inner = cram::r#async::io::writer::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository)
                    .build_with_writer(inner);
                Writer::Cram(inner)
            }
            (Format::Cram, Some(CompressionMethod::Bgzf)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CRAM cannot be compressed with BGZF",
                ));
            }
        };

        Ok(writer)
    }
}
