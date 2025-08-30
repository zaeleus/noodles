use std::path::Path;

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use tokio::{
    fs::File,
    io::{self, AsyncWrite, BufWriter},
};

use super::Writer;
use crate::variant::io::{CompressionMethod, Format};

/// An async variant writer builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
}

impl Builder {
    /// Sets the compression method of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{r#async::io::writer::Builder, io::CompressionMethod};
    /// let _builder = Builder::default().set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{r#async::io::writer::Builder, io::Format};
    /// let _builder = Builder::default().set_format(Format::Vcf);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Builds a variant writer from a path.
    ///
    /// If the format or compression method is not set, it is detected from the path extension.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant::r#async::io::writer::Builder;
    /// let _writer = Builder::default().build_from_path("out.vcf.gz").await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_path<P>(mut self, src: P) -> io::Result<Writer<File>>
    where
        P: AsRef<Path>,
    {
        use crate::variant::io::writer::builder::{
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
            .map(|file| self.build_from_writer(file))
    }

    /// Builds a variant writer from a writer.
    ///
    /// If the format is not set, a default format is used. If the compression method is not set, a
    /// default one is determined by the format.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant::r#async::io::writer::Builder;
    /// use tokio::io;
    /// let _writer = Builder::default().build_from_writer(io::sink());
    /// # Ok(())
    /// # }
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer<W>
    where
        W: AsyncWrite + Unpin,
    {
        use super::Inner;

        let format = self.format.unwrap_or(Format::Vcf);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => match format {
                Format::Vcf => None,
                Format::Bcf => Some(CompressionMethod::Bgzf),
            },
        };

        let inner = match (format, compression_method) {
            (Format::Bcf, None) => {
                Inner::BcfRaw(bcf::r#async::io::Writer::from(BufWriter::new(writer)))
            }
            (Format::Bcf, Some(CompressionMethod::Bgzf)) => {
                Inner::Bcf(bcf::r#async::io::Writer::new(writer))
            }
            (Format::Vcf, None) => {
                Inner::Vcf(vcf::r#async::io::Writer::new(BufWriter::new(writer)))
            }
            (Format::Vcf, Some(CompressionMethod::Bgzf)) => Inner::VcfGz(
                vcf::r#async::io::Writer::new(bgzf::r#async::io::Writer::new(writer)),
            ),
        };

        Writer(inner)
    }
}
