use std::{
    convert::TryFrom,
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf::{self as bgzf, index::Chunk};

use crate::index::{
    self,
    header::Format,
    reference_sequence::{self, Bin, Metadata},
    ReferenceSequence,
};

use super::{Index, MAGIC_NUMBER};

const NUL: u8 = b'\x00';

/// A tabix reader.
///
/// Consider using [`crate::read`] to read the entire index at once.
pub struct Reader<R> {
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a tabix reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_tabix as tabix;;
    /// let reader = File::open("sample.vcf.gz.tbi").map(tabix::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
        }
    }

    /// Reads the tabix index.
    ///
    /// The position of the stream is expected to be at the beginning.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_tabix as tabix;;
    /// let mut reader = File::open("sample.vcf.gz.tbi").map(tabix::Reader::new)?;
    /// let index = reader.read_index()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        read_magic(&mut self.inner)?;

        let n_ref = self.inner.read_i32::<LittleEndian>()?;

        let header = read_header(&mut self.inner)?;

        let names = read_names(&mut self.inner)?;
        let references = read_references(&mut self.inner, n_ref as usize)?;
        let n_no_coors = self.inner.read_u64::<LittleEndian>().ok();

        let mut builder = Index::builder()
            .set_header(header)
            .set_reference_sequence_names(names)
            .set_reference_sequences(references);

        if let Some(unmapped_read_count) = n_no_coors {
            builder = builder.set_unmapped_read_count(unmapped_read_count);
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
            "invalid tabix header",
        ))
    }
}

fn read_header<R>(reader: &mut R) -> io::Result<index::Header>
where
    R: Read,
{
    let format = reader.read_i32::<LittleEndian>().and_then(|n| {
        Format::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let col_seq = reader.read_i32::<LittleEndian>().and_then(|i| {
        usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let col_beg = reader.read_i32::<LittleEndian>().and_then(|i| {
        usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let col_end = reader.read_i32::<LittleEndian>().and_then(|i| {
        if i == 0 {
            Ok(None)
        } else {
            usize::try_from(i)
                .map(Some)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        }
    })?;

    let meta = reader
        .read_i32::<LittleEndian>()
        .and_then(|b| u8::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let skip = reader.read_i32::<LittleEndian>()?;

    Ok(index::Header::builder()
        .set_format(format)
        .set_reference_sequence_name_index(col_seq)
        .set_start_position_index(col_beg)
        .set_end_position_index(col_end)
        .set_line_comment_prefix(meta as u8)
        .set_line_skip_count(skip as u32)
        .build())
}

fn read_names<R>(reader: &mut R) -> io::Result<Vec<String>>
where
    R: Read,
{
    let l_nm = reader.read_i32::<LittleEndian>()?;

    let mut names = vec![0; l_nm as usize];
    reader.read_exact(&mut names)?;

    parse_names(&names)
}

fn parse_names(buf: &[u8]) -> io::Result<Vec<String>> {
    let mut names = Vec::new();
    let mut start = 0;

    loop {
        let buf = &buf[start..];

        match buf.iter().position(|&b| b == NUL) {
            Some(end) => {
                let raw_name = &buf[..end];
                let name = str::from_utf8(raw_name)
                    .map(|s| s.into())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                names.push(name);

                start += end + 1;
            }
            None => break,
        }
    }

    Ok(names)
}

fn read_references<R>(reader: &mut R, len: usize) -> io::Result<Vec<ReferenceSequence>>
where
    R: Read,
{
    let mut references = Vec::with_capacity(len);

    for _ in 0..len {
        let (bins, metadata) = read_bins(reader)?;
        let intervals = read_intervals(reader)?;
        references.push(ReferenceSequence::new(bins, intervals, metadata));
    }

    Ok(references)
}

fn read_bins<R>(reader: &mut R) -> io::Result<(Vec<Bin>, Option<Metadata>)>
where
    R: Read,
{
    use reference_sequence::bin::METADATA_ID;

    let n_bin = reader.read_i32::<LittleEndian>()?;

    let mut bins = Vec::with_capacity(n_bin as usize);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32::<LittleEndian>()?;

        if id == METADATA_ID {
            metadata = read_metadata(reader).map(Some)?;
        } else {
            let chunks = read_chunks(reader)?;
            let bin = Bin::new(id, chunks);
            bins.push(bin);
        }
    }

    Ok((bins, metadata))
}

fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: Read,
{
    let n_chunk = reader.read_i32::<LittleEndian>()?;
    let mut chunks = Vec::with_capacity(n_chunk as usize);

    for _ in 0..n_chunk {
        let cnk_beg = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        let cnk_end = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        chunks.push(Chunk::new(cnk_beg, cnk_end));
    }

    Ok(chunks)
}

fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: Read,
{
    let n_intv = reader.read_i32::<LittleEndian>()?;
    let mut intervals = Vec::with_capacity(n_intv as usize);

    for _ in 0..n_intv {
        let ioff = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioff);
    }

    Ok(intervals)
}

fn read_metadata<R>(reader: &mut R) -> io::Result<Metadata>
where
    R: Read,
{
    use reference_sequence::bin::METADATA_CHUNK_COUNT;

    let n_chunk = reader.read_u32::<LittleEndian>()?;

    if n_chunk != METADATA_CHUNK_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid metadata pseudo-bin chunk count: expected {}, got {}",
                METADATA_CHUNK_COUNT, n_chunk
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic_with_invalid_magic_number() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"TBI";
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
    fn test_read_header() -> io::Result<()> {
        let data = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1
            0x04, 0x00, 0x00, 0x00, // col_beg = 4
            0x05, 0x00, 0x00, 0x00, // col_end = 5
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
        ];

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = index::header::Builder::gff().build();
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_invalid_reference_sequence_name_index() {
        let data = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0xff, 0xff, 0xff, 0xff, // col_seq = -1
            0x04, 0x00, 0x00, 0x00, // col_beg = 4
            0x05, 0x00, 0x00, 0x00, // col_end = 5
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
        ];

        let mut reader = &data[..];
        assert!(read_header(&mut reader).is_err());
    }

    #[test]
    fn test_read_header_with_invalid_start_position_index() {
        let data = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1
            0xff, 0xff, 0xff, 0xff, // col_beg = -1
            0x05, 0x00, 0x00, 0x00, // col_end = 5
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
        ];

        let mut reader = &data[..];
        assert!(read_header(&mut reader).is_err());
    }

    #[test]
    fn test_read_header_with_invalid_end_position_index() {
        let data = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1
            0x04, 0x00, 0x00, 0x00, // col_beg = 4
            0xff, 0xff, 0xff, 0xff, // col_end = -1
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
        ];

        let mut reader = &data[..];
        assert!(read_header(&mut reader).is_err());
    }

    #[test]
    fn test_read_header_with_invalid_line_comment_prefix() {
        let data = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1
            0x04, 0x00, 0x00, 0x00, // col_beg = 4
            0x05, 0x00, 0x00, 0x00, // col_end = 5
            0x5c, 0xf3, 0x01, 0x00, // meta = '🍜'
            0x00, 0x00, 0x00, 0x00, // skip = 0
        ];

        let mut reader = &data[..];
        assert!(read_header(&mut reader).is_err());
    }

    #[test]
    fn test_parse_names() -> io::Result<()> {
        let data = b"noodles\x00tabix\x00";
        let actual = parse_names(&data[..])?;
        let expected = vec![String::from("noodles"), String::from("tabix")];
        assert_eq!(actual, expected);

        let data = b"";
        assert!(parse_names(&data[..])?.is_empty());

        Ok(())
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
}
