pub mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

use std::{fs::File, io, path::Path};

use bit_vec::BitVec;
use noodles_bgzf::VirtualPosition;

use self::index::reference::{bin::Chunk, Bin};

pub static MAGIC_NUMBER: &[u8] = b"BAI\x01";

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
        let bin_index = bin.bin() as usize;

        if bin_index < region_bins.len() && region_bins[bin_index] {
            query_bins.push(bin);
        }
    }

    query_bins
}

/// Merges a list of chunks into a list of non-overlapping chunks.
///
/// This is the same as calling [`optimize_chunks`] with a `min_offset` of 0.
pub fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    optimize_chunks(chunks, VirtualPosition::default())
}

/// Optimizes a list of chunks into a list of non-overlapping chunks.
///
/// Unlike [`merge_chunks`], `min_offset` (typically from the linear index) is given to remove
/// chunks that cannot be in the query.
pub fn optimize_chunks(chunks: &[Chunk], min_offset: VirtualPosition) -> Vec<Chunk> {
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

    fn build_chunks() -> Vec<Chunk> {
        vec![
            Chunk::new(VirtualPosition::from(2), VirtualPosition::from(5)),
            Chunk::new(VirtualPosition::from(3), VirtualPosition::from(4)),
            Chunk::new(VirtualPosition::from(5), VirtualPosition::from(7)),
            Chunk::new(VirtualPosition::from(9), VirtualPosition::from(12)),
            Chunk::new(VirtualPosition::from(10), VirtualPosition::from(15)),
            Chunk::new(VirtualPosition::from(16), VirtualPosition::from(21)),
        ]
    }

    #[test]
    fn test_merge_chunks() {
        let chunks = build_chunks();
        let actual = merge_chunks(&chunks);

        let expected = [
            Chunk::new(VirtualPosition::from(2), VirtualPosition::from(7)),
            Chunk::new(VirtualPosition::from(9), VirtualPosition::from(15)),
            Chunk::new(VirtualPosition::from(16), VirtualPosition::from(21)),
        ];

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
        let chunks = build_chunks();
        let actual = optimize_chunks(&chunks, VirtualPosition::from(10));

        let expected = [
            Chunk::new(VirtualPosition::from(9), VirtualPosition::from(15)),
            Chunk::new(VirtualPosition::from(16), VirtualPosition::from(21)),
        ];

        assert_eq!(actual, expected);
    }
}
