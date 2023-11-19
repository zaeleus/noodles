use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf as bgzf;
use noodles_csi::{
    index::{
        header::ReferenceSequenceNames,
        reference_sequence::{bin::Chunk, Bin, Metadata},
        Header, ReferenceSequence,
    },
    BinningIndex,
};

use super::{Index, MAGIC_NUMBER};

/// A tabix writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: bgzf::Writer<W>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a tabix writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let writer = tabix::Writer::new(Vec::new());
    /// ```
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::Writer::new(writer),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let writer = tabix::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        self.inner.get_ref()
    }

    /// Attempts to finish the output stream.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_tabix as tabix;
    /// let mut writer = tabix::Writer::new(Vec::new());
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
        self.inner.try_finish()
    }

    /// Writes a tabix index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// use noodles_tabix as tabix;
    ///
    /// let mut writer = tabix::Writer::new(Vec::new());
    ///
    /// let index = csi::Index::builder()
    ///     .set_header(csi::index::Header::default())
    ///     .build();
    ///
    /// writer.write_index(&index)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_index(&mut self.inner, index)
    }
}

fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: Write,
{
    write_magic(writer)?;

    let n_ref = i32::try_from(index.reference_sequences().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_ref)?;

    let header = index
        .header()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing tabix header"))?;
    write_header(writer, header)?;

    write_reference_sequences(writer, index.reference_sequences())?;

    if let Some(n_no_coor) = index.unplaced_unmapped_record_count() {
        writer.write_u64::<LittleEndian>(n_no_coor)?;
    }

    Ok(())
}

fn write_magic<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(MAGIC_NUMBER)
}

fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: Write,
{
    let format = i32::from(header.format());
    writer.write_i32::<LittleEndian>(format)?;

    let reference_sequence_name_index = header
        .reference_sequence_name_index()
        .checked_add(1)
        .expect("attempt to add with overflow");
    let col_seq = i32::try_from(reference_sequence_name_index)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(col_seq)?;

    let start_position_index = header
        .start_position_index()
        .checked_add(1)
        .expect("attempt to add with overflow");
    let col_beg = i32::try_from(start_position_index)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(col_beg)?;

    let col_end = header.end_position_index().map_or(Ok(0), |mut i| {
        i = i.checked_add(1).expect("attempt to add with overflow");
        i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;
    writer.write_i32::<LittleEndian>(col_end)?;

    let meta = i32::from(header.line_comment_prefix());
    writer.write_i32::<LittleEndian>(meta)?;

    let skip = i32::try_from(header.line_skip_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(skip)?;

    write_reference_sequence_names(writer, header.reference_sequence_names())?;

    Ok(())
}

fn write_reference_sequence_names<W>(
    writer: &mut W,
    reference_sequence_names: &ReferenceSequenceNames,
) -> io::Result<()>
where
    W: Write,
{
    const NUL: u8 = 0x00;

    // Add 1 for each trailing NUL.
    let len = reference_sequence_names
        .iter()
        .map(|n| n.len() + 1)
        .sum::<usize>();
    let l_nm = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(l_nm)?;

    for reference_sequence_name in reference_sequence_names {
        writer.write_all(reference_sequence_name.as_bytes())?;
        writer.write_u8(NUL)?;
    }

    Ok(())
}

fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &[ReferenceSequence],
) -> io::Result<()>
where
    W: Write,
{
    for reference_sequence in reference_sequences {
        write_reference_sequence(writer, reference_sequence)?;
    }

    Ok(())
}

pub fn write_reference_sequence<W>(writer: &mut W, reference: &ReferenceSequence) -> io::Result<()>
where
    W: Write,
{
    let mut n_bin = i32::try_from(reference.bins().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if reference.metadata().is_some() {
        n_bin = n_bin
            .checked_add(1)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "n_bin overflow"))?;
    }

    writer.write_i32::<LittleEndian>(n_bin)?;

    for (&id, bin) in reference.bins() {
        write_bin(writer, id, bin)?;
    }

    if let Some(metadata) = reference.metadata() {
        write_metadata(writer, metadata)?;
    }

    let n_intv = i32::try_from(reference.linear_index().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_intv)?;

    for interval in reference.linear_index() {
        let ioff = u64::from(*interval);
        writer.write_u64::<LittleEndian>(ioff)?;
    }

    Ok(())
}

pub fn write_bin<W>(writer: &mut W, id: usize, bin: &Bin) -> io::Result<()>
where
    W: Write,
{
    let id = u32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(id)?;
    write_chunks(writer, bin.chunks())?;
    Ok(())
}

fn write_chunks<W>(writer: &mut W, chunks: &[Chunk]) -> io::Result<()>
where
    W: Write,
{
    let n_chunk =
        i32::try_from(chunks.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_chunk)?;

    for chunk in chunks {
        write_chunk(writer, chunk)?;
    }

    Ok(())
}

fn write_chunk<W>(writer: &mut W, chunk: &Chunk) -> io::Result<()>
where
    W: Write,
{
    let cnk_beg = u64::from(chunk.start());
    writer.write_u64::<LittleEndian>(cnk_beg)?;

    let cnk_end = u64::from(chunk.end());
    writer.write_u64::<LittleEndian>(cnk_end)?;

    Ok(())
}

fn write_metadata<W>(writer: &mut W, metadata: &Metadata) -> io::Result<()>
where
    W: Write,
{
    use super::index::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);
    const METADATA_CHUNK_COUNT: usize = 2;

    let bin_id =
        u32::try_from(METADATA_ID).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(bin_id)?;

    let n_chunk = u32::try_from(METADATA_CHUNK_COUNT)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n_chunk)?;

    let ref_beg = u64::from(metadata.start_position());
    writer.write_u64::<LittleEndian>(ref_beg)?;

    let ref_end = u64::from(metadata.end_position());
    writer.write_u64::<LittleEndian>(ref_end)?;

    let n_mapped = metadata.mapped_record_count();
    writer.write_u64::<LittleEndian>(n_mapped)?;

    let n_unmapped = metadata.unmapped_record_count();
    writer.write_u64::<LittleEndian>(n_unmapped)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;
    use noodles_csi as csi;

    use super::*;

    #[test]
    fn test_write_index() -> io::Result<()> {
        let chunks = vec![Chunk::new(
            bgzf::VirtualPosition::from(509268599425),
            bgzf::VirtualPosition::from(509268599570),
        )];
        let bins = [(16385, Bin::new(chunks))].into_iter().collect();
        let intervals = vec![bgzf::VirtualPosition::from(337)];
        let references = vec![ReferenceSequence::new(bins, intervals, None)];

        let reference_sequence_names = [String::from("sq0"), String::from("sq1")]
            .into_iter()
            .collect();

        let header = csi::index::Header::builder()
            .set_reference_sequence_names(reference_sequence_names)
            .build();

        let index = Index::builder()
            .set_header(header)
            .set_reference_sequences(references)
            .build();

        let mut buf = Vec::new();
        write_index(&mut buf, &index)?;

        let mut expected = Vec::new();
        // magic
        expected.write_all(MAGIC_NUMBER)?;
        // n_ref
        expected.write_i32::<LittleEndian>(1)?;
        // format
        expected.write_i32::<LittleEndian>(0)?;
        // col_seq
        expected.write_i32::<LittleEndian>(1)?;
        // col_beg
        expected.write_i32::<LittleEndian>(4)?;
        // col_end
        expected.write_i32::<LittleEndian>(5)?;
        // meta
        expected.write_i32::<LittleEndian>(i32::from(b'#'))?;
        // skip
        expected.write_i32::<LittleEndian>(0)?;
        // l_nm
        expected.write_i32::<LittleEndian>(8)?;
        // names
        expected.write_all(b"sq0\x00sq1\x00")?;
        // n_bin
        expected.write_u32::<LittleEndian>(1)?;
        // bin
        expected.write_u32::<LittleEndian>(16385)?;
        // n_chunk
        expected.write_u32::<LittleEndian>(1)?;
        // chunk_beg
        expected.write_u64::<LittleEndian>(509268599425)?;
        // chunk_end
        expected.write_u64::<LittleEndian>(509268599570)?;
        // n_intv
        expected.write_u32::<LittleEndian>(1)?;
        // ioffset
        expected.write_u64::<LittleEndian>(337)?;

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_metadata() -> io::Result<()> {
        let metadata = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        let mut buf = Vec::new();
        write_metadata(&mut buf, &metadata)?;

        let expected = [
            0x4a, 0x92, 0x00, 0x00, // bin = 37450
            0x02, 0x00, 0x00, 0x00, // chunks = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_unmapped = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
