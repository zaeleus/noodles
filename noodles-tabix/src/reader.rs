use std::{
    convert::TryFrom,
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

use crate::index::{
    self,
    header::Format,
    reference_sequence::{bin::Chunk, Bin},
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

    let meta = reader.read_i32::<LittleEndian>()?;
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
        let bins = read_bins(reader)?;
        let intervals = read_intervals(reader)?;
        references.push(ReferenceSequence::new(bins, intervals));
    }

    Ok(references)
}

fn read_bins<R>(reader: &mut R) -> io::Result<Vec<Bin>>
where
    R: Read,
{
    let n_bin = reader.read_i32::<LittleEndian>()?;
    let mut bins = Vec::with_capacity(n_bin as usize);

    for _ in 0..n_bin {
        let bin = reader.read_u32::<LittleEndian>()?;
        let chunks = read_chunks(reader)?;
        bins.push(Bin::new(bin, chunks));
    }

    Ok(bins)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_index_with_invalid_magic_number() {
        let data = [];
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_index().is_err());

        let data = b"TBI";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_index().is_err());

        let data = b"MThd";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_index().is_err());
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
    fn test_parse_names() -> io::Result<()> {
        let data = b"noodles\x00tabix\x00";
        let actual = parse_names(&data[..])?;
        let expected = vec![String::from("noodles"), String::from("tabix")];
        assert_eq!(actual, expected);

        let data = b"";
        assert!(parse_names(&data[..])?.is_empty());

        Ok(())
    }
}
