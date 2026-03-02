//! Async CRAM reader.

mod builder;
mod container;
mod crc_reader;
pub mod header;
mod num;
mod query;
mod records;

use futures::Stream;
use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_sam as sam;
use tokio::io::{self, AsyncRead, AsyncSeek, AsyncSeekExt, SeekFrom};

pub use self::builder::Builder;
use self::{container::read_container, crc_reader::CrcReader};
use crate::{FileDefinition, crai, file_definition::Version, io::reader::Container};

/// An async CRAM reader.
pub struct Reader<R> {
    inner: R,
    reference_sequence_repository: fasta::Repository,
    version: Version,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let reader = cram::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
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
    /// use tokio::io;
    /// let mut reader = cram::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
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
    /// use tokio::io;
    /// let reader = cram::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
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
    /// use tokio::io;
    /// let reader = cram::r#async::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Builder::default().build_from_reader(inner)
    }

    /// Returns an async CRAM header reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::{fs::File, io::AsyncReadExt};
    ///
    /// let mut reader = File::open("sample.cram")
    ///     .await
    ///     .map(cram::r#async::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number().await?;
    /// header_reader.read_format_version().await?;
    /// header_reader.read_file_id().await?;
    ///
    /// let mut container_reader = header_reader.container_reader().await?;
    ///
    /// let _raw_header = {
    ///     let mut raw_sam_header_reader = container_reader.raw_sam_header_reader().await?;
    ///     let mut raw_header = String::new();
    ///     raw_sam_header_reader.read_to_string(&mut raw_header).await?;
    ///     raw_sam_header_reader.discard_to_end().await?;
    ///     raw_header
    /// };
    ///
    /// container_reader.discard_to_end().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
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
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::fs::File;
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    /// let file_definition = reader.read_file_definition().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        let file_definition = header::read_file_definition(&mut self.inner).await?;
        self.version = file_definition.version();
        self.version.validate()?;
        Ok(file_definition)
    }

    /// Reads the SAM header.
    ///
    /// The position of the stream is expected to be at the CRAM header container, i.e., directly
    /// after the file definition.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
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
    pub async fn read_file_header(&mut self) -> io::Result<sam::Header> {
        header::read_file_header(&mut self.inner, self.version).await
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
        self.read_file_header().await
    }

    /// Reads a container.
    ///
    /// This returns `None` if the container header is the EOF container header, which signals the
    /// end of the stream.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram::{self as cram, io::reader::Container};
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.cram").await.map(cram::r#async::io::Reader::new)?;
    /// reader.read_header().await?;
    ///
    /// let mut container = Container::default();
    ///
    /// while reader.read_container(&mut container).await? != 0 {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_container(&mut self, container: &mut Container) -> io::Result<usize> {
        read_container(&mut self.inner, container, self.version).await
    }

    /// Returns an (async) stream over records starting from the current (input) stream position.
    ///
    /// The (input) stream position is expected to be at the start of a container.
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
    ) -> impl Stream<Item = io::Result<sam::alignment::RecordBuf>> + 'r {
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
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram as cram;
    /// use tokio::io::{self, SeekFrom};
    /// let mut reader = cram::r#async::io::Reader::new(io::empty());
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
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use tokio::io;
    /// use noodles_cram as cram;
    /// let mut reader = cram::r#async::io::Reader::new(io::empty());
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
    ) -> io::Result<impl Stream<Item = io::Result<sam::alignment::RecordBuf>> + use<'r, 'h, 'i, R>>
    {
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
