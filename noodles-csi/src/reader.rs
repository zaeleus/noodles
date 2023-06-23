use std::{
    collections::HashMap,
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

use super::{
    index::{
        header::ReferenceSequenceNames,
        reference_sequence::{bin::Chunk, Bin, Metadata},
        Header, ReferenceSequence,
    },
    Index, MAGIC_NUMBER,
};

/// A CSI reader.
pub struct Reader<R> {
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a CSI reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_csi as csi;
    /// let reader = File::open("sample.bcf.csi").map(csi::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: bgzf::Reader::new(inner),
        }
    }

    /// Reads a CSI index.
    ///
    /// The position of the stream is expected to be at the beginning.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_csi as csi;
    /// let mut reader = File::open("sample.bcf.csi").map(csi::Reader::new)?;
    /// let index = reader.read_index();
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        read_magic(&mut self.inner)?;

        let min_shift = self.inner.read_i32::<LittleEndian>().and_then(|n| {
            u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let depth = self.inner.read_i32::<LittleEndian>().and_then(|n| {
            u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let header = read_aux(&mut self.inner)?;
        let reference_sequences = read_reference_sequences(&mut self.inner, depth)?;
        let n_no_coor = read_unplaced_unmapped_record_count(&mut self.inner)?;

        let mut builder = Index::builder()
            .set_min_shift(min_shift)
            .set_depth(depth)
            .set_reference_sequences(reference_sequences);

        if let Some(hdr) = header {
            builder = builder.set_header(hdr);
        }

        if let Some(n_no_coor) = n_no_coor {
            builder = builder.set_unplaced_unmapped_record_count(n_no_coor);
        }

        Ok(builder.build())
    }
}

fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CSI file format",
        ))
    }
}

fn read_aux<R>(reader: &mut R) -> io::Result<Option<Header>>
where
    R: Read,
{
    let l_aux = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if l_aux > 0 {
        let mut aux = vec![0; l_aux];
        reader.read_exact(&mut aux)?;

        let mut rdr = &aux[..];
        read_header(&mut rdr).map(Some)
    } else {
        Ok(None)
    }
}

pub(crate) fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: Read,
{
    use super::index::header::Format;

    let format = reader.read_i32::<LittleEndian>().and_then(|n| {
        Format::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let col_seq = reader.read_i32::<LittleEndian>().and_then(|i| {
        usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let col_beg = reader.read_i32::<LittleEndian>().and_then(|i| {
        usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let col_end = reader.read_i32::<LittleEndian>().and_then(|i| match i {
        0 => Ok(None),
        _ => usize::try_from(i)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    })?;

    let meta = reader
        .read_i32::<LittleEndian>()
        .and_then(|b| u8::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let skip = reader.read_i32::<LittleEndian>().and_then(|n| {
        u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let names = read_names(reader)?;

    Ok(Header::builder()
        .set_format(format)
        .set_reference_sequence_name_index(col_seq)
        .set_start_position_index(col_beg)
        .set_end_position_index(col_end)
        .set_line_comment_prefix(meta)
        .set_line_skip_count(skip)
        .set_reference_sequence_names(names)
        .build())
}

fn read_names<R>(reader: &mut R) -> io::Result<ReferenceSequenceNames>
where
    R: Read,
{
    let l_nm = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut names = vec![0; l_nm];
    reader.read_exact(&mut names)?;

    parse_names(&names)
}

pub(crate) fn parse_names(mut src: &[u8]) -> io::Result<ReferenceSequenceNames> {
    const NUL: u8 = 0x00;

    let mut names = ReferenceSequenceNames::new();

    while let Some(i) = src.iter().position(|&b| b == NUL) {
        let (raw_name, rest) = src.split_at(i);

        let name =
            str::from_utf8(raw_name).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        if !names.insert(name.into()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate reference sequence name: {name}"),
            ));
        }

        src = &rest[1..];
    }

    if src.is_empty() {
        Ok(names)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid reference sequence names",
        ))
    }
}

fn read_reference_sequences<R>(reader: &mut R, depth: u8) -> io::Result<Vec<ReferenceSequence>>
where
    R: Read,
{
    let n_ref = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    (0..n_ref)
        .map(|_| read_reference_sequence(reader, depth))
        .collect()
}

fn read_reference_sequence<R>(reader: &mut R, depth: u8) -> io::Result<ReferenceSequence>
where
    R: Read,
{
    let (bins, metadata) = read_bins(reader, depth)?;
    Ok(ReferenceSequence::new(bins, Vec::new(), metadata))
}

fn read_bins<R>(reader: &mut R, depth: u8) -> io::Result<(HashMap<usize, Bin>, Option<Metadata>)>
where
    R: Read,
{
    let n_bin = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = HashMap::with_capacity(n_bin);

    let metadata_id = Bin::metadata_id(depth);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32::<LittleEndian>().and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let loffset = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        if id == metadata_id {
            metadata = read_metadata(reader).map(Some)?;
        } else {
            let chunks = read_chunks(reader)?;
            let bin = Bin::new(loffset, chunks);
            // TODO: Check for duplicates.
            bins.insert(id, bin);
        }
    }

    Ok((bins, metadata))
}

fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: Read,
{
    let n_chunk = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut chunks = Vec::with_capacity(n_chunk);

    for _ in 0..n_chunk {
        let chunk_beg = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        let chunk_end = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        chunks.push(Chunk::new(chunk_beg, chunk_end));
    }

    Ok(chunks)
}

fn read_metadata<R>(reader: &mut R) -> io::Result<Metadata>
where
    R: Read,
{
    use crate::index::reference_sequence::bin::METADATA_CHUNK_COUNT;

    let n_chunk = reader.read_u32::<LittleEndian>()?;

    if n_chunk != METADATA_CHUNK_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid metadata pseudo-bin chunk count: expected {METADATA_CHUNK_COUNT}, got {n_chunk}"
            ),
        ));
    }

    let ref_beg = reader
        .read_u64::<LittleEndian>()
        .map(bgzf::VirtualPosition::from)?;

    let ref_end = reader
        .read_u64::<LittleEndian>()
        .map(bgzf::VirtualPosition::from)?;

    let n_mapped = reader.read_u64::<LittleEndian>()?;
    let n_unmapped = reader.read_u64::<LittleEndian>()?;

    Ok(Metadata::new(ref_beg, ref_end, n_mapped, n_unmapped))
}

fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> io::Result<Option<u64>>
where
    R: Read,
{
    match reader.read_u64::<LittleEndian>() {
        Ok(n) => Ok(Some(n)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic() {
        let data = b"CSI\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).is_ok());
    }

    #[test]
    fn test_read_magic_with_invalid_magic_number() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"CSI";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_read_metadata() -> io::Result<()> {
        let data = [
            0x02, 0x00, 0x00, 0x00, // n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_unmapped = 0
        ];

        let mut reader = &data[..];
        let actual = read_metadata(&mut reader)?;

        let expected = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_unplaced_unmapped_record_count() -> io::Result<()> {
        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_unplaced_unmapped_record_count(&mut reader)?, None);

        let data = [0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_unplaced_unmapped_record_count(&mut reader)?, Some(8));

        Ok(())
    }
}
