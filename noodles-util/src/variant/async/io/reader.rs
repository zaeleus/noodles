//! Async variant reader.

mod builder;
mod inner;

use futures::Stream;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncRead};

pub use self::builder::Builder;
use self::inner::Inner;
use crate::variant::Record;

/// An async variant reader.
pub struct Reader<R>(Inner<R>)
where
    R: AsyncRead;

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Reads the VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant;
    ///
    /// let src = b"##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    /// let mut reader = variant::r#async::io::Reader::new(&src[..]).await?;
    /// let _header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<vcf::Header> {
        self.0.read_header().await
    }

    /// Reads a variant record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant;
    ///
    /// let src = b"##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    /// let mut reader = variant::r#async::io::Reader::new(&src[..]).await?;
    /// let header = reader.read_header().await?;
    ///
    /// let mut record = variant::Record::default();
    ///
    /// while reader.read_record(&mut record).await? != 0 {
    ///     // ...
    /// }
    /// # Ok::<_, std::io::Error>(())
    /// }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        self.0.read_record(record).await
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_util::variant;
    ///
    /// let src = b"##fileformat=VCFv4.5\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    /// let mut reader = variant::r#async::io::Reader::new(&src[..]).await?;
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

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
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
    pub async fn new(reader: R) -> io::Result<Self> {
        Builder::default().build_from_reader(reader).await
    }
}
