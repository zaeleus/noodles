//! Async variant reader.

mod builder;

use std::pin::Pin;

use futures::{Stream, StreamExt};
use noodles_bcf as bcf;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncBufRead};

pub use self::builder::Builder;

/// An async variant reader.
pub enum Reader<R> {
    /// BCF.
    Bcf(bcf::r#async::io::Reader<R>),
    /// VCF.
    Vcf(vcf::r#async::io::Reader<R>),
}

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
        match self {
            Self::Bcf(reader) => reader.read_header().await,
            Self::Vcf(reader) => reader.read_header().await,
        }
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
        #[allow(clippy::type_complexity)]
        let records: Pin<
            Box<dyn Stream<Item = io::Result<Box<dyn vcf::variant::Record>>>>,
        > = match self {
            Reader::Bcf(reader) => Box::pin(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            })),
            Reader::Vcf(reader) => Box::pin(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn vcf::variant::Record>)
            })),
        };

        records
    }
}
