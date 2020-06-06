//! BAM index (BAI)

pub mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

use std::{fs::File, io, path::Path};

use noodles_bgzf::VirtualPosition;

use self::index::reference::{bin::Chunk, Bin};

pub(crate) static MAGIC_NUMBER: &[u8] = b"BAI\x01";

pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_header()?;
    reader.read_index()
}

/// Merges a list of chunks into a list of non-overlapping chunks.
///
/// This is the same as calling [`optimize_chunks`] with a `min_offset` of 0.
///
/// [`optimize_chunks`]: fn.optimize_chunks.html
pub fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    optimize_chunks(chunks, VirtualPosition::default())
}

/// Optimizes a list of chunks into a list of non-overlapping chunks.
///
/// Unlike [`merge_chunks`], `min_offset` (typically from the linear index) is given to remove
/// chunks that cannot be in the query.
///
/// [`merge_chunks`]: fn.merge_chunks.html
pub fn optimize_chunks(chunks: &[Chunk], min_offset: VirtualPosition) -> Vec<Chunk> {
    let mut chunks: Vec<_> = chunks
        .iter()
        .filter(|c| c.end() > min_offset)
        .copied()
        .collect();

    if chunks.is_empty() {
        return chunks;
    }

    chunks.sort_unstable_by_key(|c| c.start());

    // At worst, no chunks are merged, and the resulting list will be the same size as the input.
    let mut merged_chunks = Vec::with_capacity(chunks.len());

    // `chunks` is guaranteed to be non-empty.
    let mut current_chunk = chunks[0];

    for next_chunk in chunks.iter().skip(1) {
        if next_chunk.start() > current_chunk.end() {
            merged_chunks.push(current_chunk);
            current_chunk = *next_chunk;
        } else if current_chunk.end() < next_chunk.end() {
            current_chunk = Chunk::new(current_chunk.start(), next_chunk.end());
        }
    }

    merged_chunks.push(current_chunk);

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
