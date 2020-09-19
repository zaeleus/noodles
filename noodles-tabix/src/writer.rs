use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use super::{
    index::{
        reference_sequence::{bin::Chunk, Bin},
        ReferenceSequence,
    },
    Index, MAGIC_NUMBER,
};

const NUL: u8 = b'\x00';

/// A tabix writer.
pub struct Writer<W> {
    inner: W,
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
    pub fn new(inner: W) -> Self {
        Self { inner }
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
        &self.inner
    }

    /// Writes a tabix index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_tabix as tabix;
    /// let index = tabix::Index::default();
    /// let mut writer = tabix::Writer::new(Vec::new());
    /// writer.write_index(&index)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_magic(&mut self.inner)?;

        let n_ref = index.reference_sequences().len() as i32;
        self.inner.write_i32::<LittleEndian>(n_ref)?;

        let format = i32::from(index.format());
        self.inner.write_i32::<LittleEndian>(format)?;

        let col_seq = index.reference_sequence_name_index() as i32;
        self.inner.write_i32::<LittleEndian>(col_seq)?;

        let col_beg = index.start_position_index() as i32;
        self.inner.write_i32::<LittleEndian>(col_beg)?;

        let col_end = index.end_position_index() as i32;
        self.inner.write_i32::<LittleEndian>(col_end)?;

        let meta = i32::from(index.line_comment_prefix());
        self.inner.write_i32::<LittleEndian>(meta)?;

        let skip = index.header_line_count() as i32;
        self.inner.write_i32::<LittleEndian>(skip)?;

        // Add 1 for each trailing nul.
        let l_nm = index
            .reference_sequence_names()
            .iter()
            .map(|n| n.len() + 1)
            .sum::<usize>() as i32;
        self.inner.write_i32::<LittleEndian>(l_nm)?;

        for reference_sequence_name in index.reference_sequence_names() {
            self.inner.write_all(reference_sequence_name.as_bytes())?;
            self.inner.write_u8(NUL)?;
        }

        for reference_sequence in index.reference_sequences() {
            write_reference_sequence(&mut self.inner, reference_sequence)?;
        }

        if let Some(n_no_coor) = index.unmapped_read_count() {
            self.inner.write_u64::<LittleEndian>(n_no_coor)?;
        }

        Ok(())
    }
}

pub fn write_magic<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(MAGIC_NUMBER)
}

pub fn write_reference_sequence<W>(writer: &mut W, reference: &ReferenceSequence) -> io::Result<()>
where
    W: Write,
{
    let n_bin = reference.bins().len() as i32;
    writer.write_i32::<LittleEndian>(n_bin)?;

    for bin in reference.bins() {
        write_bin(writer, bin)?;
    }

    let n_intv = reference.intervals().len() as i32;
    writer.write_i32::<LittleEndian>(n_intv)?;

    for interval in reference.intervals() {
        let ioff = u64::from(*interval);
        writer.write_u64::<LittleEndian>(ioff)?;
    }

    Ok(())
}

pub fn write_bin<W>(writer: &mut W, bin: &Bin) -> io::Result<()>
where
    W: Write,
{
    writer.write_u32::<LittleEndian>(bin.id())?;

    let n_chunk = bin.chunks().len() as i32;
    writer.write_i32::<LittleEndian>(n_chunk)?;

    for chunk in bin.chunks() {
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

#[cfg(test)]
mod tests {
    use std::io::BufWriter;

    use noodles_bgzf as bgzf;

    use crate::index::Format;

    use super::*;

    #[test]
    fn test_write_index() -> io::Result<()> {
        let chunks = vec![Chunk::new(
            bgzf::VirtualPosition::from(509268599425),
            bgzf::VirtualPosition::from(509268599570),
        )];
        let bins = vec![Bin::new(16385, chunks)];
        let intervals = vec![bgzf::VirtualPosition::from(337)];
        let references = vec![ReferenceSequence::new(bins, intervals)];
        let index = Index::builder()
            .set_format(Format::Vcf)
            .set_reference_sequence_name_index(1)
            .set_start_position_index(4)
            .set_end_position_index(5)
            .set_line_comment_prefix(b'#')
            .set_header_line_count(0)
            .set_reference_sequence_names(vec![String::from("sq0"), String::from("sq1")])
            .set_reference_sequences(references)
            .build();

        let mut actual_writer = Writer::new(Vec::new());
        actual_writer.write_index(&index)?;

        let mut expected_writer = BufWriter::new(Vec::new());
        // magic
        expected_writer.write_all(MAGIC_NUMBER)?;
        // n_ref
        expected_writer.write_i32::<LittleEndian>(1)?;
        // format
        expected_writer.write_i32::<LittleEndian>(2)?;
        // col_seq
        expected_writer.write_i32::<LittleEndian>(1)?;
        // col_beg
        expected_writer.write_i32::<LittleEndian>(4)?;
        // col_end
        expected_writer.write_i32::<LittleEndian>(5)?;
        // meta
        expected_writer.write_i32::<LittleEndian>(i32::from(b'#'))?;
        // skip
        expected_writer.write_i32::<LittleEndian>(0)?;
        // l_nm
        expected_writer.write_i32::<LittleEndian>(8)?;
        // names
        expected_writer.write_all(b"sq0\x00sq1\x00")?;
        // n_bin
        expected_writer.write_u32::<LittleEndian>(1)?;
        // bin
        expected_writer.write_u32::<LittleEndian>(16385)?;
        // n_chunk
        expected_writer.write_u32::<LittleEndian>(1)?;
        // chunk_beg
        expected_writer.write_u64::<LittleEndian>(509268599425)?;
        // chunk_end
        expected_writer.write_u64::<LittleEndian>(509268599570)?;
        // n_intv
        expected_writer.write_u32::<LittleEndian>(1)?;
        // ioffset
        expected_writer.write_u64::<LittleEndian>(337)?;
        expected_writer.flush()?;

        let actual = actual_writer.get_ref();
        let expected = expected_writer.get_ref();

        assert_eq!(actual, expected);

        Ok(())
    }
}
