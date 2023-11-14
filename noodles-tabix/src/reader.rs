use std::{
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_csi::index::{
    header::{Format, ReferenceSequenceNames},
    reference_sequence::{bin::Chunk, Bin, Metadata},
    Header, ReferenceSequence,
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

        let n_ref = self.inner.read_i32::<LittleEndian>().and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let header = read_header(&mut self.inner)?;

        let references = read_references(&mut self.inner, n_ref)?;
        let n_no_coor = read_unplaced_unmapped_record_count(&mut self.inner)?;

        let mut builder = Index::builder()
            .set_header(header)
            .set_reference_sequences(references);

        if let Some(unplaced_unmapped_record_count) = n_no_coor {
            builder = builder.set_unplaced_unmapped_record_count(unplaced_unmapped_record_count);
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

fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: Read,
{
    let format = reader.read_i32::<LittleEndian>().and_then(|n| {
        Format::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let col_seq = reader.read_i32::<LittleEndian>().and_then(|i| {
        usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|n| {
                if i > 0 {
                    Ok(n - 1)
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid col_seq",
                    ))
                }
            })
    })?;

    let col_beg = reader.read_i32::<LittleEndian>().and_then(|i| {
        usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|n| {
                if i > 0 {
                    Ok(n - 1)
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid col_bed",
                    ))
                }
            })
    })?;

    let col_end = reader.read_i32::<LittleEndian>().and_then(|i| {
        if i == 0 {
            Ok(None)
        } else {
            usize::try_from(i)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                .and_then(|n| {
                    if n > 0 {
                        Ok(n - 1)
                    } else {
                        Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "invalid col_end",
                        ))
                    }
                })
                .map(Some)
        }
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

pub(crate) fn parse_names(mut buf: &[u8]) -> io::Result<ReferenceSequenceNames> {
    let mut names = ReferenceSequenceNames::new();

    while let Some(i) = buf.iter().position(|&b| b == NUL) {
        let (raw_name, rest) = buf.split_at(i);

        let name =
            str::from_utf8(raw_name).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        if !names.insert(name.into()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate reference sequence name: {name}"),
            ));
        }

        buf = &rest[1..];
    }

    if buf.is_empty() {
        Ok(names)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid reference sequence names",
        ))
    }
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

fn read_bins<R>(reader: &mut R) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)>
where
    R: Read,
{
    use super::index::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);

    let n_bin = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = IndexMap::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32::<LittleEndian>().and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let is_duplicate = if id == METADATA_ID {
            let m = read_metadata(reader)?;
            metadata.replace(m).is_some()
        } else {
            let chunks = read_chunks(reader)?;
            let bin = Bin::new(bgzf::VirtualPosition::default(), chunks);
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

fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: Read,
{
    let n_chunk = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut chunks = Vec::with_capacity(n_chunk);

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
    let n_intv = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

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
    const METADATA_CHUNK_COUNT: usize = 2;

    let n_chunk = reader.read_u32::<LittleEndian>().and_then(|n| {
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
        use noodles_csi::index::header;

        let data = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1
            0x04, 0x00, 0x00, 0x00, // col_beg = 4
            0x05, 0x00, 0x00, 0x00, // col_end = 5
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
            0x00, 0x00, 0x00, 0x00, // l_nm = 0
        ];

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = header::Builder::gff().build();
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
            0x00, 0x00, 0x00, 0x00, // l_nm = 0
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
            0x00, 0x00, 0x00, 0x00, // l_nm = 0
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
            0x00, 0x00, 0x00, 0x00, // l_nm = 0
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
            0x00, 0x00, 0x00, 0x00, // l_nm = 0
        ];

        let mut reader = &data[..];
        assert!(read_header(&mut reader).is_err());
    }

    #[test]
    fn test_parse_names() -> io::Result<()> {
        let data = b"noodles\x00tabix\x00";
        let actual = parse_names(&data[..])?;
        let expected: ReferenceSequenceNames = [String::from("noodles"), String::from("tabix")]
            .into_iter()
            .collect();
        assert_eq!(actual, expected);

        let data = b"";
        assert!(parse_names(&data[..])?.is_empty());

        Ok(())
    }

    #[test]
    fn test_parse_names_with_duplicate_name() {
        let data = b"sq0\x00sq0\x00";

        assert!(matches!(
            parse_names(data),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }

    #[test]
    fn test_parse_names_with_trailing_data() {
        let data = b"sq0\x00sq1\x00sq2";

        assert!(matches!(
            parse_names(data),
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
