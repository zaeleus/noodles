use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;

use crate::index::{
    reference::{bin::Chunk, Bin},
    Reference,
};

use super::{Index, MAGIC_NUMBER};

pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    pub fn read_index(&mut self) -> io::Result<Index> {
        read_magic(&mut self.inner)?;

        let n_ref = self.inner.read_i32::<LittleEndian>()?;
        let format = self.inner.read_i32::<LittleEndian>()?;
        let col_seq = self.inner.read_i32::<LittleEndian>()?;
        let col_beg = self.inner.read_i32::<LittleEndian>()?;
        let col_end = self.inner.read_i32::<LittleEndian>()?;
        let meta = self.inner.read_i32::<LittleEndian>()?;
        let skip = self.inner.read_i32::<LittleEndian>()?;
        let names = read_names(&mut self.inner)?;
        let references = read_references(&mut self.inner, n_ref as usize)?;
        let n_no_coors = self.inner.read_u64::<LittleEndian>().ok();

        Ok(Index::new(
            format,
            col_seq as usize,
            col_beg as usize,
            col_end as usize,
            meta,
            skip as u32,
            names,
            references,
            n_no_coors,
        ))
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

fn read_names<R>(reader: &mut R) -> io::Result<Vec<String>>
where
    R: Read,
{
    let l_nm = reader.read_i32::<LittleEndian>()?;

    let mut names = vec![0; l_nm as usize];
    reader.read_exact(&mut names)?;

    let list = Vec::new();

    Ok(list)
}

fn read_references<R>(reader: &mut R, len: usize) -> io::Result<Vec<Reference>>
where
    R: Read,
{
    let mut references = Vec::with_capacity(len);

    for _ in 0..len {
        let bins = read_bins(reader)?;
        let intervals = read_intervals(reader)?;
        references.push(Reference::new(bins, intervals));
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
}
