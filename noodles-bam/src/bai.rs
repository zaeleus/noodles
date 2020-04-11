mod reader;

pub use self::reader::Reader;

use std::{fs::File, io, path::Path};

use bit_vec::BitVec;

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

const WINDOW_SIZE: u64 = 16384;

#[derive(Debug)]
pub struct Reference {
    pub bins: Vec<Bin>,
    pub intervals: Vec<Interval>,
}

impl Reference {
    pub fn min_offset(&self, start: u64) -> u64 {
        let i = (start / WINDOW_SIZE) as usize;
        self.intervals[i].ioffset
    }
}

#[derive(Debug)]
pub struct Index {
    pub references: Vec<Reference>,
    pub n_no_coor: Option<u64>,
}

pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_header()?;
    reader.read_index()
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

const MIN_OFFSET: u64 = 0;

/// Merges a list of chunks into a list of non-overlapping chunks.
///
/// This is the same as calling [`optimize_chunks`] with a `min_offset` of 0.
pub fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    optimize_chunks(chunks, MIN_OFFSET)
}

/// Optimizes a list of chunks into a list of non-overlapping chunks.
///
/// Unlike [`merge_chunks`], `min_offset` (typically from the linear index) is given to remove
/// chunks that cannot be in the query.
pub fn optimize_chunks(chunks: &[Chunk], min_offset: u64) -> Vec<Chunk> {
    let mut chunks: Vec<_> = chunks
        .iter()
        .filter(|c| c.end() > min_offset)
        .cloned()
        .collect();

    if chunks.is_empty() {
        return chunks;
    }

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

    #[test]
    fn test_optimize_chunks() {
        let chunks = [
            Chunk::new(2, 5),
            Chunk::new(3, 4),
            Chunk::new(5, 7),
            Chunk::new(9, 12),
            Chunk::new(10, 15),
            Chunk::new(16, 21),
        ];

        let actual = optimize_chunks(&chunks, 10);
        let expected = [Chunk::new(9, 15), Chunk::new(16, 21)];

        assert_eq!(actual, expected);
    }
}
