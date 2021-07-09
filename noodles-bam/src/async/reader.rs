use std::convert::TryFrom;

use noodles_bgzf as bgzf;
use noodles_sam::header::{ReferenceSequence, ReferenceSequences};
use pin_project_lite::pin_project;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{reader::bytes_with_nul_to_string, Record, MAGIC_NUMBER};

pin_project! {
    /// An async BAM reader.
    pub struct Reader<R>
    where
        R: AsyncRead,
    {
        #[pin]
        inner: bgzf::AsyncReader<R>,
    }
}

impl<R> Reader<R>
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
        Self {
            inner: bgzf::AsyncReader::new(reader),
        }
    }

    /// Reads the raw SAM header.
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
    /// The position of the stream is expected to be directly after the header.
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
    /// The stream is expected to be directly after the reference sequences or at the start of
    /// another record.
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
    /// reader.read_reference_sequences().await?;
    ///
    /// let mut record = bam::Record::default();
    /// reader.read_record(&mut record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record).await
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
        reference_sequences.insert(reference_sequence.name().into(), reference_sequence);
    }

    Ok(reference_sequences)
}

async fn read_reference_sequence<R>(reader: &mut R) -> io::Result<ReferenceSequence>
where
    R: AsyncRead + Unpin,
{
    let l_name = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut c_name = vec![0; l_name];
    reader.read_exact(&mut c_name).await?;

    let name = bytes_with_nul_to_string(&c_name)?;
    let l_ref = reader.read_u32_le().await.and_then(|len| {
        i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    ReferenceSequence::new(name, l_ref).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

async fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let block_size = match reader.read_u32_le().await {
        Ok(bs) => usize::try_from(bs).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    };

    record.resize(block_size);
    reader.read_exact(record).await?;

    Ok(block_size)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_magic() {
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
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            dbg!(read_magic(&mut reader).await),
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
        let data = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequences(&mut reader).await?;

        let expected: ReferenceSequences = vec![("sq0", 8)]
            .into_iter()
            .map(|(name, len)| ReferenceSequence::new(name, len).map(|rs| (name.into(), rs)))
            .collect::<Result<_, _>>()?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_record() -> io::Result<()> {
        let data = [
            0x04, 0x00, 0x00, 0x00, // block_size = 4
            0x6e, 0x64, 0x6c, 0x73, // ...
        ];

        let mut reader = &data[..];
        let mut record = Record::default();
        let block_size = read_record(&mut reader, &mut record).await?;

        assert_eq!(block_size, 4);
        assert_eq!(&record[..], &data[4..]);

        Ok(())
    }
}
