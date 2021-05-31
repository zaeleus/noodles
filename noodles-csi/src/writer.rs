use std::{
    convert::TryFrom,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_bgzf::{self as bgzf, index::Chunk};

use super::{
    index::{reference_sequence::Bin, ReferenceSequence},
    Index, MAGIC_NUMBER,
};

/// A CSI writer.
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
    /// Creates a CSI writer.
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::Writer::new(writer),
        }
    }

    /// Writes a coordinate-sorted index (CSI).
    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_magic(&mut self.inner)?;

        let min_shift = index.min_shift();
        self.inner.write_i32::<LittleEndian>(min_shift)?;

        let depth = index.depth();
        self.inner.write_i32::<LittleEndian>(depth)?;

        write_aux(&mut self.inner, index.aux())?;
        write_reference_sequences(&mut self.inner, index.reference_sequences())?;

        if let Some(n_no_coor) = index.unmapped_read_count() {
            self.inner.write_u64::<LittleEndian>(n_no_coor)?;
        }

        Ok(())
    }
}

fn write_magic<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(MAGIC_NUMBER)
}

fn write_aux<W>(writer: &mut W, aux: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let l_aux =
        i32::try_from(aux.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(l_aux)?;

    writer.write_all(aux)?;

    Ok(())
}

fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &[ReferenceSequence],
) -> io::Result<()>
where
    W: Write,
{
    let n_ref = i32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_ref)?;

    for reference_sequence in reference_sequences {
        write_bins(writer, reference_sequence.bins())?;
    }

    Ok(())
}

fn write_bins<W>(writer: &mut W, bins: &[Bin]) -> io::Result<()>
where
    W: Write,
{
    let n_bin =
        i32::try_from(bins.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_bin)?;

    for bin in bins {
        let bin_id = bin.id();
        writer.write_u32::<LittleEndian>(bin_id)?;

        write_chunks(writer, bin.chunks())?;
    }

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
        let chunk_beg = u64::from(chunk.start());
        writer.write_u64::<LittleEndian>(chunk_beg)?;

        let chunk_end = u64::from(chunk.start());
        writer.write_u64::<LittleEndian>(chunk_end)?;
    }

    Ok(())
}
