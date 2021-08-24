mod container;
mod data_container;
mod num;

use tokio::io::{self, AsyncRead, AsyncReadExt, AsyncSeek, AsyncSeekExt, SeekFrom};

use crate::{file_definition::Version, DataContainer, FileDefinition};

/// An async CRAM reader.
pub struct Reader<R> {
    inner: R,
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
    /// let reader = cram::AsyncReader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
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
    /// let mut reader = File::open("sample.cram").await.map(cram::AsyncReader::new)?;
    /// let file_definition = reader.read_file_definition().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        read_magic_number(&mut self.inner).await?;

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
    /// let mut reader = File::open("sample.cram").await.map(cram::AsyncReader::new)?;
    /// reader.read_file_definition().await?;
    ///
    /// let header = reader.read_file_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_file_header(&mut self) -> io::Result<String> {
        use self::container::read_container;
        use crate::reader::read_file_header_block;

        let container = read_container(&mut self.inner).await?;

        if let Some(block) = container.blocks().first() {
            read_file_header_block(block)
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid header container: missing block for SAM header",
            ))
        }
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
    /// let mut reader = File::open("sample.cram").await.map(cram::AsyncReader::new)?;
    /// reader.read_file_definition().await?;
    /// reader.read_file_header().await?;
    ///
    /// while let Some(container) = reader.read_data_container().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_data_container(&mut self) -> io::Result<Option<DataContainer>> {
        use self::data_container::read_data_container;

        read_data_container(&mut self.inner).await
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
    /// let mut reader = cram::AsyncReader::new(Cursor::new(Vec::new()));
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
    /// let mut reader = cram::AsyncReader::new(Cursor::new(Vec::new()));
    /// assert_eq!(reader.position().await?, 0);
    /// # Ok(())
    /// # }
    /// ```
    pub async fn position(&mut self) -> io::Result<u64> {
        self.inner.seek(SeekFrom::Current(0)).await
    }
}

async fn read_magic_number<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    use crate::MAGIC_NUMBER;

    let mut magic = [0; 4];
    reader.read_exact(&mut magic).await?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CRAM header",
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
    async fn test_read_magic_number() {
        let data = b"CRAM";
        let mut reader = &data[..];
        assert!(read_magic_number(&mut reader).await.is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"BAM\x01";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic_number(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

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
