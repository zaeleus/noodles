use std::io::{self, Read};

use bit_vec::BitVec;
use byteorder::{LittleEndian, ReadBytesExt};

pub static MAGIC_NUMBER: &[u8] = b"BAI\x01";

#[derive(Debug)]
pub struct Interval {
    ioffset: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Chunk {
    chunk_beg: u64,
    chunk_end: u64,
}

impl Chunk {
    pub fn new(start: u64, end: u64) -> Self {
        Self {
            chunk_beg: start,
            chunk_end: end,
        }
    }

    pub fn start(&self) -> u64 {
        self.chunk_beg
    }

    pub fn start_mut(&mut self) -> &mut u64 {
        &mut self.chunk_beg
    }

    pub fn end(&self) -> u64 {
        self.chunk_end
    }

    pub fn end_mut(&mut self) -> &mut u64 {
        &mut self.chunk_end
    }
}

#[derive(Debug)]
pub struct Bin {
    bin: u32,
    chunks: Vec<Chunk>,
}

impl Bin {
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }
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
    pub fn new(inner: R) -> Self {
        Self { inner }
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

const MAX_BINS: usize = ((1 << 18) - 1) / 7 + 1;

// See § 5.3 in SAMv1.pdf (accessed 2019-11-15).
pub fn region_to_bins(start: usize, end: usize) -> BitVec {
    let ranges = [
        (1 + (start >> 26), 1 + (end >> 26)),
        (9 + (start >> 23), 9 + (end >> 23)),
        (73 + (start >> 20), 73 + (end >> 20)),
        (585 + (start >> 17), 585 + (end >> 17)),
        (4681 + (start >> 14), 4681 + (end >> 14)),
    ];

    let mut bins = BitVec::from_elem(MAX_BINS, false);

    bins.set(0, true);

    for (start, end) in &ranges {
        for k in *start..=*end {
            bins.set(k, true);
        }
    }

    bins
}

pub fn query(bins: &[Bin], start: u64, end: u64) -> Vec<&Bin> {
    let region_bins = region_to_bins(start as usize, end as usize);

    let mut query_bins = Vec::new();

    for bin in bins {
        let bin_index = bin.bin as usize;

        if bin_index < region_bins.len() && region_bins[bin_index] {
            query_bins.push(bin);
        }
    }

    query_bins
}

pub fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    if chunks.is_empty() {
        return Vec::new();
    }

    let mut chunks = chunks.to_vec();
    chunks.sort_unstable_by_key(|c| c.start());

    let mut merged_chunks = Vec::with_capacity(chunks.len());
    merged_chunks.push(chunks[0].clone());

    for b in chunks {
        let a = merged_chunks.last_mut().expect("list cannot be empty");

        if b.start() > a.end() {
            merged_chunks.push(b);
            continue;
        }

        if a.end() < b.end() {
            *a.end_mut() = b.end();
        }
    }

    merged_chunks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_chunks() {
        let chunks = [
            Chunk::new(2, 5),
            Chunk::new(3, 4),
            Chunk::new(5, 7),
            Chunk::new(9, 12),
            Chunk::new(10, 15),
            Chunk::new(16, 21),
        ];

        let actual = merge_chunks(&chunks);
        let expected = [Chunk::new(2, 7), Chunk::new(9, 15), Chunk::new(16, 21)];

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_merge_chunks_with_empty_list() {
        let chunks = Vec::new();
        let merged_chunks = merge_chunks(&chunks);
        assert!(merged_chunks.is_empty());
    }
}
