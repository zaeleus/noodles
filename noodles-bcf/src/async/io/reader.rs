mod header;
mod query;
mod record;

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncRead, AsyncReadExt, AsyncSeek};

use self::{header::read_header, query::query, record::read_record};
use crate::{
    header::{string_maps::ContigStringMap, StringMaps},
    Record,
};

/// An async BCF reader.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// #
/// # #[tokio::main]
/// # async fn main() -> io::Result<()> {
/// use futures::TryStreamExt;
/// use noodles_bcf as bcf;
/// use tokio::fs::File;
///
/// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
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
pub struct Reader<R> {
    inner: R,
    string_maps: StringMaps,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let data = [];
    /// let reader = bcf::r#async::io::Reader::from(&data[..]);
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
    /// use noodles_bcf as bcf;
    /// let data = [];
    /// let mut reader = bcf::r#async::io::Reader::from(&data[..]);
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
    /// use noodles_bcf as bcf;
    /// let data = [];
    /// let reader = bcf::r#async::io::Reader::from(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Returns the string maps.
    ///
    /// This is only built after reading the header using [`Self::read_header`].
    pub fn string_maps(&self) -> &StringMaps {
        &self.string_maps
    }

    /// Reads the VCF header.
    ///
    /// The BCF magic number is checked, and the file format version is discarded.
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
    /// use noodles_bcf as bcf;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<vcf::Header> {
        read_magic(&mut self.inner).await?;
        read_format_version(&mut self.inner).await?;
        let (header, string_maps) = read_header(&mut self.inner).await?;
        self.string_maps = string_maps;
        Ok(header)
    }

    /// Reads a single record without decoding (most of) its feilds.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergonomic to read records using a stream (see [`Self::records`]), but using
    /// this method directly allows the reuse of a single [`Record`] buffer.
    ///
    /// If successful, the record size is returned. If a record size of 0 is returned, the stream
    /// reached EOF.
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bcf as bcf;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
    /// reader.read_header().await?;
    ///
    /// let mut record = bcf::Record::default();
    /// reader.read_record(&mut record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record).await
    }

    /// Returns an (async) stream over lazy records starting from the current (input) stream
    /// position.
    ///
    /// The (input) stream is expected to be directly after the header or at the start of another
    /// record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_bcf as bcf;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
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
    pub fn records(&mut self) -> impl Stream<Item = io::Result<Record>> + '_ {
        Box::pin(stream::try_unfold(
            (&mut self.inner, Record::default()),
            |(mut reader, mut record)| async {
                read_record(&mut reader, &mut record)
                    .await
                    .map(|n| match n {
                        0 => None,
                        _ => Some((record.clone(), (reader, record))),
                    })
            },
        ))
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async BCF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let data = [];
    /// let reader = bcf::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self::from(bgzf::AsyncReader::new(inner))
    }

    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// use noodles_bgzf as bgzf;
    ///
    /// let data = Vec::new();
    /// let reader = bcf::r#async::io::Reader::new(&data[..]);
    /// let virtual_position = reader.virtual_position();
    ///
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn virtual_position(&self) -> bgzf::VirtualPosition {
        self.inner.virtual_position()
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Seeks the underlying BGZF reader to the given virtual position.
    ///
    /// Virtual positions typically come from an associated index file.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bcf as bcf;
    /// use noodles_bgzf as bgzf;
    ///
    /// let data = [];
    /// let mut reader = bcf::r#async::io::Reader::new(Cursor::new(data));
    ///
    /// let virtual_position = bgzf::VirtualPosition::default();
    /// reader.seek(virtual_position).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek(pos).await
    }

    /// Returns a stream over records that intersect the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_bcf as bcf;
    /// use noodles_core::Region;
    /// use noodles_csi as csi;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf").await.map(bcf::r#async::io::Reader::new)?;
    /// reader.read_header().await?;
    ///
    /// let string_maps = reader.string_maps().clone();
    ///
    /// let index = csi::r#async::read("sample.bcf.csi").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(string_maps.contigs(), &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<I>(
        &mut self,
        contig_string_map: &ContigStringMap,
        index: &I,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + '_>
    where
        I: BinningIndex,
    {
        use crate::io::reader::resolve_region;

        let reference_sequence_id = resolve_region(contig_string_map, region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(query(
            self,
            chunks,
            reference_sequence_id,
            region.interval(),
        ))
    }
}

impl<R> From<R> for Reader<R> {
    fn from(inner: R) -> Self {
        Self {
            inner,
            string_maps: StringMaps::default(),
        }
    }
}

async fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    use crate::MAGIC_NUMBER;

    let mut magic = [0; 3];
    reader.read_exact(&mut magic).await?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BCF header",
        ))
    }
}

async fn read_format_version<R>(reader: &mut R) -> io::Result<(u8, u8)>
where
    R: AsyncRead + Unpin,
{
    let major_version = reader.read_u8().await?;
    let minor_version = reader.read_u8().await?;

    Ok((major_version, minor_version))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_magic() {
        let data = b"BCF";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).await.is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"BAM";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[tokio::test]
    async fn test_read_format_version() -> io::Result<()> {
        let data = [0x02, 0x01];
        let mut reader = &data[..];
        assert_eq!(read_format_version(&mut reader).await?, (2, 1));
        Ok(())
    }
}
