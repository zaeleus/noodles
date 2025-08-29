use std::path::Path;

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use tokio::{
    fs::File,
    io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, BufReader},
};

use super::Reader;
use crate::variant::io::{CompressionMethod, Format};

/// An async variant reader builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
}

impl Builder {
    /// Sets the compression method of the input.
    ///
    /// By default, the compression method is autodetected on build. This can be used to override
    /// it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{r#async::io::reader::Builder, io::CompressionMethod};
    /// let _builder = Builder::default().set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression);
        self
    }

    /// Sets the format of the input.
    ///
    /// By default, the format is autodetected on build. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{r#async::io::reader::Builder, io::Format};
    /// let _builder = Builder::default().set_format(Format::Vcf);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Builds an async variant reader from a path.
    ///
    /// By default, the format and compression method will be autodetected. This can be overridden
    /// by using [`Self::set_format`] and [`Self::set_compression_method`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant::r#async::io::reader::Builder;
    /// let _reader = Builder::default().build_from_path("samples.vcf").await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_path<P>(
        self,
        src: P,
    ) -> io::Result<Reader<Box<dyn AsyncBufRead + Unpin>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(src).await?;
        self.build_from_reader(file).await
    }

    /// Builds an async variant reader from a reader.
    ///
    /// By default, the format and compression method will be autodetected. This can be overridden
    /// by using [`Self::set_format`] and [`Self::set_compression_method`].
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant::r#async::io::reader::Builder;
    /// use tokio::io;
    /// let reader = Builder::default().build_from_reader(io::empty()).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_reader<'a, R>(
        self,
        reader: R,
    ) -> io::Result<Reader<Box<dyn AsyncBufRead + Unpin + 'a>>>
    where
        R: AsyncRead + Unpin + 'a,
    {
        use super::Inner;
        use crate::variant::io::reader::builder::{detect_compression_method, detect_format};

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
            (Format::Vcf, None) => {
                let inner: Box<dyn AsyncBufRead + Unpin> = Box::new(reader);
                Inner::Vcf(vcf::r#async::io::Reader::new(inner))
            }
            (Format::Vcf, Some(CompressionMethod::Bgzf)) => {
                let decoder: Box<dyn AsyncBufRead + Unpin> =
                    Box::new(bgzf::r#async::io::Reader::new(reader));
                Inner::Vcf(vcf::r#async::io::Reader::new(decoder))
            }
            (Format::Bcf, None) => {
                let inner: Box<dyn AsyncBufRead + Unpin> = Box::new(reader);
                Inner::Bcf(bcf::r#async::io::Reader::from(inner))
            }
            (Format::Bcf, Some(CompressionMethod::Bgzf)) => {
                let decoder: Box<dyn AsyncBufRead + Unpin> =
                    Box::new(bgzf::r#async::io::Reader::new(reader));
                Inner::Bcf(bcf::r#async::io::Reader::from(decoder))
            }
        };

        Ok(Reader(inner))
    }
}
