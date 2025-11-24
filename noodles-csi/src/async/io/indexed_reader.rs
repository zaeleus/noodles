use futures::Stream;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::{Query, filtered_indexed_records};
use crate::{BinningIndex, io::indexed_records::Record};

/// An async indexed reader.
pub struct IndexedReader<R, I> {
    inner: R,
    index: I,
}

impl<R, I> IndexedReader<bgzf::r#async::io::Reader<R>, I>
where
    R: AsyncRead,
    I: BinningIndex,
{
    /// Creates an async indexed reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    ///
    /// let index = csi::Index::default();
    /// let reader = csi::r#async::io::IndexedReader::new(io::empty(), index);
    /// ```
    pub fn new(inner: R, index: I) -> Self {
        Self {
            inner: bgzf::r#async::io::Reader::new(inner),
            index,
        }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    ///
    /// let index = csi::Index::default();
    /// let reader = csi::r#async::io::IndexedReader::new(io::empty(), index);
    /// let inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::r#async::io::Reader<R> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    ///
    /// let index = csi::Index::default();
    /// let mut reader = csi::r#async::io::IndexedReader::new(io::empty(), index);
    /// let inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::r#async::io::Reader<R> {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    ///
    /// let index = csi::Index::default();
    /// let reader = csi::r#async::io::IndexedReader::new(io::empty(), index);
    /// let inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::r#async::io::Reader<R> {
        self.inner
    }

    /// Returns the associated index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    ///
    /// let index = csi::Index::default();
    /// let reader = csi::r#async::io::IndexedReader::new(io::empty(), index);
    /// let index = reader.index();
    /// ```
    pub fn index(&self) -> &I {
        &self.index
    }
}

impl<R, I> IndexedReader<bgzf::r#async::io::Reader<R>, I>
where
    R: AsyncRead + AsyncSeek + Unpin,
    I: BinningIndex,
{
    /// Returns a stream over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_csi as csi;
    /// use tokio::fs::File;
    ///
    /// let index = csi::r#async::fs::read("sample.bed.gz.csi").await?;
    ///
    /// let mut reader = File::open("sample.bed.gz")
    ///     .await
    ///     .map(|f| csi::r#async::io::IndexedReader::new(f, index))?;
    ///
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query(&mut self, region: &Region) -> io::Result<impl Stream<Item = io::Result<Record>>> {
        let header = self
            .index
            .header()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing index header"))?;

        let reference_sequence_id = header
            .reference_sequence_names()
            .get_index_of(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "missing reference sequence name",
                )
            })?;

        let chunks = self.index.query(reference_sequence_id, region.interval())?;

        Ok(filtered_indexed_records(
            Query::new(&mut self.inner, chunks),
            header,
            region,
        ))
    }
}
