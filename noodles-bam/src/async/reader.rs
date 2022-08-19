mod query;
mod record;

use bytes::BytesMut;
use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_sam::{
    alignment::Record,
    header::{
        record::value::{map::ReferenceSequence, Map},
        ReferenceSequences,
    },
};
use tokio::io::{self, AsyncRead, AsyncReadExt, AsyncSeek};

use self::{query::query, record::read_record};
use crate::{
    lazy,
    reader::{bytes_with_nul_to_string, resolve_region},
    MAGIC_NUMBER,
};

/// An async BAM reader.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// #
/// # #[tokio::main]
/// # async fn main() -> io::Result<()> {
/// use futures::TryStreamExt;
/// use noodles_bam as bam;
/// use tokio::fs::File;
///
/// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
/// reader.read_header().await?;
/// reader.read_reference_sequences().await?;
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
    buf: BytesMut,
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
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::AsyncReader::from(&data[..]);
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
    /// use noodles_bam as bam;
    /// let data = [];
    /// let mut reader = bam::AsyncReader::from(&data[..]);
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
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::AsyncReader::from(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads the raw SAM header.
    ///
    /// The BAM magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
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
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
    /// let header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<String> {
        read_magic(&mut self.inner).await?;
        read_header(&mut self.inner).await
    }

    /// Reads the binary reference sequences after the SAM header.
    ///
    /// This is not the same as the `@SQ` records in the SAM header. A BAM has a list of reference
    /// sequences containing name and length tuples after the SAM header and before the list of
    /// records.
    ///
    /// The position of the stream is expected to be directly after the header.
    ///
    /// This returns a reference sequence dictionary ([`noodles_sam::header::ReferenceSequences`]),
    /// which can be used to build a minimal [`noodles_sam::Header`] if the SAM header is empty.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
    /// reader.read_header().await?;
    /// let reference_sequences = reader.read_reference_sequences().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_reference_sequences(&mut self) -> io::Result<ReferenceSequences> {
        read_reference_sequences(&mut self.inner).await
    }

    /// Reads a single record.
    ///
    /// The record block size (`bs`) is read from the underlying stream, and `bs` bytes are read
    /// into the given record buffer.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// It is more ergonomic to read records using a stream (see [`Self::records`] and
    /// [`Self::query`]), but using this method directly allows the reuse of a single [`Record`]
    /// buffer.
    ///
    /// If successful, the record block size is returned. If a block size of 0 is returned, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_sam::alignment::Record;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
    /// reader.read_header().await?;
    /// reader.read_reference_sequences().await?;
    ///
    /// let mut record = Record::default();
    /// reader.read_record(&mut record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, &mut self.buf, record).await
    }

    /// Reads a single record without eagerly decoding its fields.
    ///
    /// The record block size (`bs`) is read from the underlying (input) stream and `bs` bytes are
    /// read into the lazy record's buffer. No fields are decoded, meaning the record is not
    /// necessarily valid. However, the structure of the byte stream is guaranteed to be
    /// record-like.
    ///
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
    ///
    /// If successful, the record block size is returned. If a block size of 0 is returned, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_sam::alignment::Record;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
    /// reader.read_header().await?;
    /// reader.read_reference_sequences().await?;
    ///
    /// let mut record = bam::lazy::Record::default();
    /// reader.read_lazy_record(&mut record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_lazy_record(&mut self, record: &mut lazy::Record) -> io::Result<usize> {
        use self::record::read_block_size;

        let block_size = match read_block_size(&mut self.inner).await? {
            0 => return Ok(0),
            n => n,
        };

        record.buf.resize(block_size, 0);
        self.inner.read_exact(&mut record.buf).await?;

        record.index()?;

        Ok(block_size)
    }

    /// Returns an (async) stream over records starting from the current (input) stream position.
    ///
    /// The (input) stream is expected to be directly after the reference sequences or at the start
    /// of another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
    /// reader.read_header().await?;
    /// reader.read_reference_sequences().await?;
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
            (&mut self.inner, &mut self.buf, Record::default()),
            move |(reader, buf, mut record)| async move {
                read_record(reader, buf, &mut record)
                    .await
                    .map(|n| match n {
                        0 => None,
                        _ => Some((record.clone(), (reader, buf, record))),
                    })
            },
        ))
    }

    /// Returns a stream over lazy records.
    ///
    /// The (input) stream is expected to be directly after the reference sequences or at the start
    /// of another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_bam as bam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
    /// reader.read_header().await?;
    /// reader.read_reference_sequences().await?;
    ///
    /// let mut records = reader.lazy_records();
    ///
    /// while let Some(record) = records.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn lazy_records(&mut self) -> impl Stream<Item = io::Result<lazy::Record>> + '_ {
        Box::pin(stream::try_unfold(
            (self, lazy::Record::default()),
            |(this, mut record)| async {
                this.read_lazy_record(&mut record).await.map(|n| match n {
                    0 => None,
                    _ => Some((record.clone(), (this, record))),
                })
            },
        ))
    }
}

impl<R> Reader<bgzf::AsyncReader<R>>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async BAM reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::AsyncReader::new(&data[..]);
    /// ```
    pub fn new(reader: R) -> Self {
        Self::from(bgzf::AsyncReader::new(reader))
    }

    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use noodles_bgzf as bgzf;
    ///
    /// let data = Vec::new();
    /// let reader = bam::AsyncReader::new(&data[..]);
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
    /// Virtual positions typically come from the associated BAM index file.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_bgzf as bgzf;
    ///
    /// let data = [];
    /// let mut reader = bam::AsyncReader::new(Cursor::new(data));
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
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>>{
    /// use futures::TryStreamExt;
    /// use noodles_bam::{self as bam, bai};
    /// use noodles_core::Region;
    /// use noodles_sam as sam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::AsyncReader::new)?;
    /// let header: sam::Header = reader.read_header().await?.parse()?;
    ///
    /// let reference_sequences = header.reference_sequences();
    /// let index = bai::r#async::read("sample.bam.bai").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(reference_sequences, &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<I>(
        &mut self,
        reference_sequences: &ReferenceSequences,
        index: &I,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Record>> + '_>
    where
        I: BinningIndex,
    {
        let reference_sequence_id = resolve_region(reference_sequences, region)?;

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
            buf: BytesMut::new(),
        }
    }
}

async fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic).await?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM header",
        ))
    }
}

async fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncRead + Unpin,
{
    let l_text = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut text = vec![0; l_text];
    reader.read_exact(&mut text).await?;

    // § 4.2 The BAM format (2021-06-03): "Plain header text in SAM; not necessarily
    // NUL-terminated".
    bytes_with_nul_to_string(&text).or_else(|_| {
        String::from_utf8(text).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

async fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
where
    R: AsyncRead + Unpin,
{
    let n_ref = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = ReferenceSequences::with_capacity(n_ref);

    for _ in 0..n_ref {
        let reference_sequence = read_reference_sequence(reader).await?;
        let name = reference_sequence.name().to_string();
        reference_sequences.insert(name, reference_sequence);
    }

    Ok(reference_sequences)
}

async fn read_reference_sequence<R>(reader: &mut R) -> io::Result<Map<ReferenceSequence>>
where
    R: AsyncRead + Unpin,
{
    let l_name = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut c_name = vec![0; l_name];
    reader.read_exact(&mut c_name).await?;

    let name = bytes_with_nul_to_string(&c_name).and_then(|name| {
        name.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let l_ref = reader.read_u32_le().await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    Map::<ReferenceSequence>::new(name, l_ref)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use noodles_sam as sam;

    use super::*;

    #[tokio::test]
    async fn test_read_magic() {
        let data = b"BAM\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).await.is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[tokio::test]
    async fn test_read_header() -> io::Result<()> {
        let expected = "@HD\tVN:1.6\n";

        let data_len = expected.len() as u32;
        let mut data = data_len.to_le_bytes().to_vec();
        data.extend(expected.as_bytes());

        let mut reader = &data[..];
        let actual = read_header(&mut reader).await?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        use sam::header::record::value::map::reference_sequence::Name;

        let data = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequences(&mut reader).await?;

        let expected: ReferenceSequences = [("sq0".parse()?, 8)]
            .into_iter()
            .map(|(name, len): (Name, usize)| {
                let sn = name.to_string();
                Map::<ReferenceSequence>::new(name, len).map(|rs| (sn, rs))
            })
            .collect::<Result<_, _>>()?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
