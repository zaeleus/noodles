//! BAM index (BAI) and fields.
//!
//! A BAM index (BAI) is used with an associated coordinate-sorted BAM file that allows random
//! access to records, e.g., [querying].
//!
//! The index contains a list of reference sequences parallel to the one defined in the BAM file.
//! Each indexed reference sequence has a calculated set of hierarchical bins at different
//! granularities. The bins then define a list of physical file positions in the BAM to search for
//! overlapping records.
//!
//! When reading entire BAM files sequentially, a BAM index is not necessary.
//!
//! [querying]: crate::Reader::query
//!
//! # Examples
//!
//! ## Reading a BAM index
//!
//! ```no_run
//! # use std::io;
//! use noodles_bam::bai;
//! let index = bai::read("sample.bam.bai")?;
//! # Ok::<(), io::Error>(())
//! ```

pub mod index;
mod reader;
mod writer;

pub use self::{index::Index, reader::Reader, writer::Writer};

use std::{fs::File, io, path::Path};

use noodles_bgzf::{index::Chunk, VirtualPosition};

use self::index::reference_sequence::Bin;

static MAGIC_NUMBER: &[u8] = b"BAI\x01";

/// Reads the entire contents of a BAM index.
///
/// This is a convenience function and is equivalent to opening the file at the given path, reading
/// the header, and reading the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_bam::bai;
/// let index = bai::read("sample.bam.bai")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_header()?;
    reader.read_index()
}

/// Writes a BAM index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path, writing
/// the header, and writing the index.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_bam::bai;
/// let index = bai::Index::default();
/// bai::write("sample.bam.bai", &index)?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &Index) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(Writer::new)?;
    writer.write_header()?;
    writer.write_index(index)
}

/// Merges a list of chunks into a list of non-overlapping chunks.
///
/// This is the same as calling [`optimize_chunks`] with a `min_offset` of 0.
///
/// # Examples
///
/// ```
/// use noodles_bam::bai;
/// use noodles_bgzf::{self as bgzf, index::Chunk};
///
/// let chunks = [
///     Chunk::new(bgzf::VirtualPosition::from(2), bgzf::VirtualPosition::from(3)),
///     Chunk::new(bgzf::VirtualPosition::from(5), bgzf::VirtualPosition::from(8)),
///     Chunk::new(bgzf::VirtualPosition::from(7), bgzf::VirtualPosition::from(13)),
///     Chunk::new(bgzf::VirtualPosition::from(21), bgzf::VirtualPosition::from(34)),
/// ];
///
/// let actual = bai::merge_chunks(&chunks);
///
/// let expected = [
///     Chunk::new(bgzf::VirtualPosition::from(2), bgzf::VirtualPosition::from(3)),
///     Chunk::new(bgzf::VirtualPosition::from(5), bgzf::VirtualPosition::from(13)),
///     Chunk::new(bgzf::VirtualPosition::from(21), bgzf::VirtualPosition::from(34)),
/// ];
///
/// assert_eq!(actual, expected);
/// ```
pub fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    optimize_chunks(chunks, VirtualPosition::default())
}

/// Optimizes a list of chunks into a list of non-overlapping chunks.
///
/// Unlike [`merge_chunks`], `min_offset` (typically from the linear index) is given to remove
/// chunks that cannot be in the query.
///
/// # Examples
///
/// ```
/// use noodles_bam::bai;
/// use noodles_bgzf::{self as bgzf, index::Chunk};
///
/// let chunks = [
///     Chunk::new(bgzf::VirtualPosition::from(2), bgzf::VirtualPosition::from(3)),
///     Chunk::new(bgzf::VirtualPosition::from(5), bgzf::VirtualPosition::from(8)),
///     Chunk::new(bgzf::VirtualPosition::from(7), bgzf::VirtualPosition::from(13)),
///     Chunk::new(bgzf::VirtualPosition::from(21), bgzf::VirtualPosition::from(34)),
/// ];
/// let min_offset = bgzf::VirtualPosition::from(5);
///
/// let actual = bai::optimize_chunks(&chunks, min_offset);
///
/// let expected = [
///     Chunk::new(bgzf::VirtualPosition::from(5), bgzf::VirtualPosition::from(13)),
///     Chunk::new(bgzf::VirtualPosition::from(21), bgzf::VirtualPosition::from(34)),
/// ];
///
/// assert_eq!(actual, expected);
/// ```
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
