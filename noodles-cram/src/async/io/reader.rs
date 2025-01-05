//! Async CRAM reader.

mod builder;
mod crc_reader;
mod data_container;
mod header;
mod header_container;
mod num;
mod query;
mod records;

use bytes::BytesMut;
use futures::Stream;
use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_sam as sam;
use tokio::io::{self, AsyncRead, AsyncReadExt, AsyncSeek, AsyncSeekExt, SeekFrom};

pub use self::builder::Builder;
use self::crc_reader::CrcReader;
use crate::{crai, file_definition::Version, DataContainer, FileDefinition, Record};

/// An async CRAM reader.
pub struct Reader<R> {
    inner: R,
    reference_sequence_repository: fasta::Repository,
    buf: BytesMut,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async CRAM reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let data = [];
    /// let reader = cram::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Builder::default().build_from_reader(inner)
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let data = [];
    /// let reader = cram::r#async::io::Reader::new(&data[..]);
    /// assert!(reader.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let data = [];
    /// let mut reader = cram::r#async::io::Reader::new(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let data = [];
    /// let reader = cram::r#async::io::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    fn reference_sequence_repository(&self) -> &fasta::Repository {
        &self.reference_sequence_repository
    }

    /// Reads the CRAM file definition.
    ///
    /// This also checks the magic number.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::fs::File;
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    /// let file_definition = reader.read_file_definition().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        use crate::io::reader::header::magic_number;

        header::read_magic_number(&mut self.inner)
            .await
            .and_then(magic_number::validate)?;

        let format = read_format(&mut self.inner).await?;
        let file_id = read_file_id(&mut self.inner).await?;

        Ok(FileDefinition::new(format, file_id))
    }

    /// Reads the raw SAM header.
    ///
    /// The position of the stream is expected to be at the CRAM header container, i.e., directly
    /// after the file definition.
    ///
    /// This returns the raw SAM header as a [`String`]. It can subsequently be parsed as a
    /// [`noodles_sam::Header`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    /// reader.read_file_definition().await?;
    ///
    /// let header = reader.read_file_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_file_header(&mut self) -> io::Result<String> {
        use self::header_container::read_header_container;
        read_header_container(&mut self.inner, &mut self.buf).await
    }

    /// Reads the SAM header.
    ///
    /// This verifies the CRAM magic number, discards the file definition, and reads and parses the
    /// file header as a SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::fs::File;
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    /// let _header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<sam::Header> {
        self.read_file_definition().await?;

        self.read_file_header().await.and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    /// Reads a data container.
    ///
    /// This returns `None` if the container header is the EOF container header, which signals the
    /// end of the stream.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    /// reader.read_header().await?;
    ///
    /// while let Some(container) = reader.read_data_container().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_data_container(&mut self) -> io::Result<Option<DataContainer>> {
        use self::data_container::read_data_container;

        read_data_container(&mut self.inner, &mut self.buf).await
    }

    /// Returns an (async) stream over records starting from the current (input) stream position.
    ///
    /// The (input) stream position is expected to be at the start of a data container.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_cram as cram;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    ///
    /// let mut records = reader.records(&header);
    ///
    /// while let Some(record) = records.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Stream<Item = io::Result<Record>> + 'r {
        use self::records::records;

        records(self, header)
    }
}

impl<R> Reader<R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Seeks the underlying reader to the given position.
    ///
    /// Positions typically come from an associated CRAM index file.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use std::io::{Cursor, SeekFrom};
    /// use noodles_cram as cram;
    /// let mut reader = cram::r#async::io::Reader::new(Cursor::new(Vec::new()));
    /// reader.seek(SeekFrom::Start(0)).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.inner.seek(pos).await
    }

    /// Returns the current position of the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use std::io::{Cursor, SeekFrom};
    /// use noodles_cram as cram;
    /// let mut reader = cram::r#async::io::Reader::new(Cursor::new(Vec::new()));
    /// assert_eq!(reader.position().await?, 0);
    /// # Ok(())
    /// # }
    /// ```
    pub async fn position(&mut self) -> io::Result<u64> {
        self.inner.stream_position().await
    }

    /// Returns a stream over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_core::Region;
    /// use noodles_cram::{self as cram, crai};
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    ///
    /// let index = crai::r#async::read("sample.cram.crai").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<'r, 'h: 'r, 'i: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &'i crai::Index,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + 'r> {
        use self::query::query;

        let reference_sequence_id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid reference sequence name",
                )
            })?;

        Ok(query(
            self,
            header,
            index,
            reference_sequence_id,
            region.interval(),
        ))
    }
}

async fn read_format<R>(reader: &mut R) -> io::Result<Version>
where
    R: AsyncRead + Unpin,
{
    let major = reader.read_u8().await?;
    let minor = reader.read_u8().await?;
    Ok(Version::new(major, minor))
}

async fn read_file_id<R>(reader: &mut R) -> io::Result<[u8; 20]>
where
    R: AsyncRead + Unpin,
{
    let mut file_id = [0; 20];
    reader.read_exact(&mut file_id).await?;
    Ok(file_id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_format() -> io::Result<()> {
        let data = [0x03, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_format(&mut reader).await?, Version::new(3, 0));
        Ok(())
    }

    #[tokio::test]
    async fn test_read_file_id() -> io::Result<()> {
        let data = [
            0x00, 0xac, 0x24, 0xf8, 0xc4, 0x2d, 0xc2, 0xa5, 0x56, 0xa0, 0x85, 0x1c, 0xa5, 0xef,
            0xf0, 0xfc, 0x6d, 0x40, 0x33, 0x4d,
        ];

        let mut reader = &data[..];
        assert_eq!(read_file_id(&mut reader).await?, data);

        Ok(())
    }
}
