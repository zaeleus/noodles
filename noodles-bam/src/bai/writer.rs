use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use super::{
    index::{
        reference::{bin::Chunk, Bin},
        Reference,
    },
    Index, MAGIC_NUMBER,
};

pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    pub fn write_header(&mut self) -> io::Result<()> {
        self.inner.write_all(MAGIC_NUMBER)
    }

    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        let n_ref = index.references().len() as u32;
        self.inner.write_u32::<LittleEndian>(n_ref)?;

        for reference in index.references() {
            write_reference(&mut self.inner, reference)?;
        }

        if let Some(n_no_coor) = index.unplaced_unmapped_read_count() {
            self.inner.write_u64::<LittleEndian>(n_no_coor)?;
        }

        Ok(())
    }
}

fn write_reference<W>(writer: &mut W, reference: &Reference) -> io::Result<()>
where
    W: Write,
{
    let n_bin = reference.bins().len() as u32;
    writer.write_u32::<LittleEndian>(n_bin)?;

    for bin in reference.bins() {
        write_bin(writer, bin)?;
    }

    let n_intv = reference.intervals().len() as u32;
    writer.write_u32::<LittleEndian>(n_intv)?;

    for interval in reference.intervals() {
        let ioffset = u64::from(*interval);
        writer.write_u64::<LittleEndian>(ioffset)?;
    }

    Ok(())
}

fn write_bin<W>(writer: &mut W, bin: &Bin) -> io::Result<()>
where
    W: Write,
{
    writer.write_u32::<LittleEndian>(bin.bin())?;

    let n_chunk = bin.chunks().len() as u32;
    writer.write_u32::<LittleEndian>(n_chunk)?;

    for chunk in bin.chunks() {
        write_chunk(writer, chunk)?;
    }

    Ok(())
}

fn write_chunk<W>(writer: &mut W, chunk: &Chunk) -> io::Result<()>
where
    W: Write,
{
    let chunk_beg = u64::from(chunk.start());
    writer.write_u64::<LittleEndian>(chunk_beg)?;

    let chunk_end = u64::from(chunk.end());
    writer.write_u64::<LittleEndian>(chunk_end)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::BufWriter;

    use noodles_bgzf as bgzf;

    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());
        writer.write_header()?;
        assert_eq!(writer.get_ref(), b"BAI\x01");
        Ok(())
    }

    #[test]
    fn test_write_index() -> io::Result<()> {
        let chunks = vec![Chunk::new(
            bgzf::VirtualPosition::from(509268599425),
            bgzf::VirtualPosition::from(509268599570),
        )];
        let bins = vec![Bin::new(16385, chunks)];
        let intervals = vec![bgzf::VirtualPosition::from(337)];
        let references = vec![Reference::new(bins, intervals)];
        let index = Index::new(references, None);

        let mut actual_writer = Writer::new(Vec::new());
        actual_writer.write_header()?;
        actual_writer.write_index(&index)?;

        let mut expected_writer = BufWriter::new(Vec::new());
        // magic
        expected_writer.write_all(MAGIC_NUMBER)?;
        // n_ref
        expected_writer.write_u32::<LittleEndian>(1)?;
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
