//! Async variant reader.

mod builder;
mod inner;

use futures::Stream;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncBufRead};

pub use self::builder::Builder;
use self::inner::Inner;

/// An async variant reader.
pub struct Reader<R>(Inner<R>);

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    /// Reads the VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant::r#async::io::reader::Builder;
    ///
    /// let data = b"##fileformat=VCFv4.5
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// ";
    ///
    /// let mut reader = Builder::default().build_from_reader(&data[..]).await?;
    /// let _header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<vcf::Header> {
        self.0.read_header().await
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_util::variant::r#async::io::reader::Builder;
    ///
    /// let data = b"##fileformat=VCFv4.5
    /// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    /// ";
    ///
    /// let mut reader = Builder::default().build_from_reader(&data[..]).await?;
    /// reader.read_header().await?;
    ///
    /// let mut records = reader.records();
    ///
    /// while let Some(record) = records.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records(
        &mut self,
    ) -> impl Stream<Item = io::Result<Box<dyn vcf::variant::Record>>> + '_ {
        self.0.records()
    }
}

impl<'a> Reader<Box<dyn AsyncBufRead + Unpin + 'a>> {
    /// Creates a variant reader.
    ///
    /// This attempts to autodetect the compression method and format of the input.
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant;
    /// use tokio::io;
    /// let reader = variant::r#async::io::Reader::new(io::empty()).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn new<R>(reader: R) -> io::Result<Self>
    where
        R: AsyncBufRead + Unpin + 'a,
    {
        Builder::default().build_from_reader(reader).await
    }
}
