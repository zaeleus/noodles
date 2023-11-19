use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_csi::{
    self as csi,
    index::{
        reference_sequence::{bin::Chunk, Bin, Metadata},
        Header, ReferenceSequence,
    },
};
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::Index;

/// An async tabix reader.
pub struct Reader<R>
where
    R: AsyncRead,
{
    inner: bgzf::AsyncReader<R>,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async tabix reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let data = [];
    /// let reader = tabix::AsyncReader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: bgzf::AsyncReader::new(inner),
        }
    }

    /// Reads the tabix index.
    ///
    /// The position of the stream is expected to be at the beginning.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_tabix as tabix;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.vcf.gz.tbi")
    ///     .await
    ///     .map(tabix::AsyncReader::new)?;
    ///
    /// let index = reader.read_index().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        read_index(&mut self.inner).await
    }
}

async fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: AsyncRead + Unpin,
{
    read_magic(reader).await?;

    let n_ref = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let header = read_header(reader).await?;
    let reference_sequences = read_reference_sequences(reader, n_ref).await?;
    let unplaced_unmapped_record_count = read_unplaced_unmapped_record_count(reader).await?;

    let mut builder = Index::builder()
        .set_header(header)
        .set_reference_sequences(reference_sequences);

    if let Some(count) = unplaced_unmapped_record_count {
        builder = builder.set_unplaced_unmapped_record_count(count);
    }

    Ok(builder.build())
}

async fn read_magic<R>(reader: &mut R) -> io::Result<()>
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
            "invalid tabix header",
        ))
    }
}

async fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: AsyncRead + Unpin,
{
    use std::mem;

    use csi::reader::index::read_header;

    const STATIC_FIELD_COUNT: usize = 6;
    const STATIC_HEADER_SIZE: usize = mem::size_of::<i32>() * STATIC_FIELD_COUNT;

    let mut buf = vec![0; STATIC_HEADER_SIZE];
    reader.read_exact(&mut buf).await?;

    let l_nm = reader.read_i32_le().await?;

    let names_len =
        usize::try_from(l_nm).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut names = vec![0; names_len];
    reader.read_exact(&mut names).await?;

    buf.extend(l_nm.to_le_bytes());
    buf.extend(names);

    let mut buf_reader = &buf[..];
    read_header(&mut buf_reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

async fn read_reference_sequences<R>(
    reader: &mut R,
    reference_sequence_count: usize,
) -> io::Result<Vec<ReferenceSequence>>
where
    R: AsyncRead + Unpin,
{
    let mut reference_sequences = Vec::with_capacity(reference_sequence_count);

    for _ in 0..reference_sequence_count {
        let reference_sequence = read_reference_sequence(reader).await?;
        reference_sequences.push(reference_sequence);
    }

    Ok(reference_sequences)
}

async fn read_reference_sequence<R>(reader: &mut R) -> io::Result<ReferenceSequence>
where
    R: AsyncRead + Unpin,
{
    let (bins, metadata) = read_bins(reader).await?;
    let intervals = read_intervals(reader).await?;
    Ok(ReferenceSequence::new(bins, intervals, metadata))
}

async fn read_bins<R>(reader: &mut R) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)>
where
    R: AsyncRead + Unpin,
{
    use crate::index::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);

    let n_bin = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = IndexMap::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32_le().await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let is_duplicate = if id == METADATA_ID {
            let m = read_metadata(reader).await?;
            metadata.replace(m).is_some()
        } else {
            let chunks = read_chunks(reader).await?;
            let bin = Bin::new(chunks);
            bins.insert(id, bin).is_some()
        };

        if is_duplicate {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate bin ID: {id}"),
            ));
        }
    }

    Ok((bins, metadata))
}

async fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: AsyncRead + Unpin,
{
    let n_chunk = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut chunks = Vec::with_capacity(n_chunk);

    for _ in 0..n_chunk {
        let chunk = read_chunk(reader).await?;
        chunks.push(chunk);
    }

    Ok(chunks)
}

async fn read_chunk<R>(reader: &mut R) -> io::Result<Chunk>
where
    R: AsyncRead + Unpin,
{
    let cnk_beg = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    let cnk_end = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    Ok(Chunk::new(cnk_beg, cnk_end))
}

async fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: AsyncRead + Unpin,
{
    let n_intv = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

    for _ in 0..n_intv {
        let ioff = reader
            .read_u64_le()
            .await
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioff);
    }

    Ok(intervals)
}

async fn read_metadata<R>(reader: &mut R) -> io::Result<Metadata>
where
    R: AsyncRead + Unpin,
{
    const METADATA_CHUNK_COUNT: usize = 2;

    let n_chunk = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if n_chunk != METADATA_CHUNK_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid metadata pseudo-bin chunk count: expected {METADATA_CHUNK_COUNT}, got {n_chunk}"
            ),
        ));
    }

    let ref_beg = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    let ref_end = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    let n_mapped = reader.read_u64_le().await?;
    let n_unmapped = reader.read_u64_le().await?;

    Ok(Metadata::new(ref_beg, ref_end, n_mapped, n_unmapped))
}

async fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> io::Result<Option<u64>>
where
    R: AsyncRead + Unpin,
{
    match reader.read_u64_le().await {
        Ok(n_no_coor) => Ok(Some(n_no_coor)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_magic() {
        let data = b"TBI\x01";
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
        let data = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1
            0x04, 0x00, 0x00, 0x00, // col_beg = 4
            0x05, 0x00, 0x00, 0x00, // col_end = 5
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
            0x08, 0x00, 0x00, 0x00, // l_nm = 8
            b's', b'q', b'0', 0x00, // names[0] = "sq0"
            b's', b'q', b'1', 0x00, // names[1] = "sq1"
        ];

        let mut reader = &data[..];
        let actual = read_header(&mut reader).await?;

        let names = [String::from("sq0"), String::from("sq1")]
            .into_iter()
            .collect();
        let expected = csi::index::header::Builder::gff()
            .set_reference_sequence_names(names)
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
