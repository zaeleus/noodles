use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

pub static MAGIC_NUMBER: &[u8] = b"BAI\x01";

#[derive(Debug)]
pub struct Interval {
    ioffset: u64,
}

#[derive(Debug)]
pub struct Chunk {
    chunk_beg: u64,
    chunk_end: u64,
}

#[derive(Debug)]
pub struct Bin {
    bin: u32,
    chunks: Vec<Chunk>,
}

#[derive(Debug)]
pub struct Reference {
    pub bins: Vec<Bin>,
    pub intervals: Vec<Interval>,
}

#[derive(Debug)]
pub struct Index {
    pub references: Vec<Reference>,
    pub n_no_coor: Option<u64>,
}

pub struct Reader<R: Read> {
    inner: R,
}

impl<R: Read> Reader<R> {
    pub fn new(inner: R) -> Reader<R> {
        Reader { inner }
    }

    pub fn header(&mut self) -> io::Result<()> {
        let mut magic = [0; 4];
        self.inner.read_exact(&mut magic)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAI header",
            ));
        }

        Ok(())
    }

    pub fn read_index(&mut self) -> io::Result<Index> {
        let n_ref = self.inner.read_i32::<LittleEndian>()?;
        let mut references = Vec::with_capacity(n_ref as usize);
        self.read_references(&mut references)?;

        let n_no_coor = self.inner.read_u64::<LittleEndian>().ok();

        Ok(Index {
            references,
            n_no_coor,
        })
    }

    pub fn read_references(&mut self, references: &mut Vec<Reference>) -> io::Result<()> {
        let n_ref = references.capacity();

        for _ in 0..n_ref {
            let n_bin = self.inner.read_i32::<LittleEndian>()?;
            let mut bins = Vec::with_capacity(n_bin as usize);
            self.read_bins(&mut bins)?;

            let n_intv = self.inner.read_i32::<LittleEndian>()?;
            let mut intervals = Vec::with_capacity(n_intv as usize);
            self.read_intervals(&mut intervals)?;

            references.push(Reference { bins, intervals });
        }

        Ok(())
    }

    pub fn read_bins(&mut self, bins: &mut Vec<Bin>) -> io::Result<()> {
        let n_bin = bins.capacity();

        for _ in 0..n_bin {
            let bin = self.inner.read_u32::<LittleEndian>()?;

            let n_chunks = self.inner.read_i32::<LittleEndian>()?;
            let mut chunks = Vec::with_capacity(n_chunks as usize);
            self.read_chunks(&mut chunks)?;

            bins.push(Bin { bin, chunks });
        }

        Ok(())
    }

    pub fn read_chunks(&mut self, chunks: &mut Vec<Chunk>) -> io::Result<()> {
        let n_chunks = chunks.capacity();

        for _ in 0..n_chunks {
            let chunk_beg = self.inner.read_u64::<LittleEndian>()?;
            let chunk_end = self.inner.read_u64::<LittleEndian>()?;
            chunks.push(Chunk {
                chunk_beg,
                chunk_end,
            });
        }

        Ok(())
    }

    pub fn read_intervals(&mut self, intervals: &mut Vec<Interval>) -> io::Result<()> {
        let n_intervals = intervals.capacity();

        for _ in 0..n_intervals {
            let ioffset = self.inner.read_u64::<LittleEndian>()?;
            intervals.push(Interval { ioffset });
        }

        Ok(())
    }
}
