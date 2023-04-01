//! Tabix index and fields.

pub mod builder;
pub mod header;
mod indexer;

pub use self::{builder::Builder, header::Header, indexer::Indexer};

use std::io;

use noodles_core::{region::Interval, Position};
use noodles_csi::{
    binning_index::optimize_chunks,
    index::{reference_sequence::bin::Chunk, ReferenceSequence},
    BinningIndex,
};

const MIN_SHIFT: u8 = 14;
pub(crate) const DEPTH: u8 = 5;

const MAX_POSITION: Position = match Position::new((1 << (MIN_SHIFT + 3 * DEPTH)) - 1) {
    Some(position) => position,
    None => panic!(),
};

/// A tabix index.
#[derive(Debug)]
pub struct Index {
    header: Header,
    reference_sequences: Vec<ReferenceSequence>,
    unplaced_unmapped_record_count: Option<u64>,
}

impl Index {
    /// Returns a builder to create an index from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let builder = tabix::Index::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns an indexer to create an index from records.
    pub fn indexer() -> Indexer {
        Indexer::default()
    }

    /// Returns the header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let header = tabix::index::Header::default();
    /// let index = tabix::Index::builder().set_header(header.clone()).build();
    /// assert_eq!(index.header(), &header);
    /// ```
    pub fn header(&self) -> &Header {
        &self.header
    }
}

impl BinningIndex for Index {
    type ReferenceSequence = ReferenceSequence;

    /// Returns a list of indexed reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::BinningIndex;
    /// use noodles_tabix as tabix;
    /// let index = tabix::Index::default();
    /// assert!(index.reference_sequences().is_empty());
    /// ```
    fn reference_sequences(&self) -> &[Self::ReferenceSequence] {
        &self.reference_sequences
    }

    /// Returns the number of unplaced, unmapped records in the associated file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::BinningIndex;
    /// use noodles_tabix as tabix;
    /// let index = tabix::Index::default();
    /// assert!(index.unplaced_unmapped_record_count().is_none());
    /// ```
    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        self.unplaced_unmapped_record_count
    }

    fn query<I>(&self, reference_sequence_id: usize, interval: I) -> io::Result<Vec<Chunk>>
    where
        I: Into<Interval>,
    {
        let reference_sequence = self
            .reference_sequences()
            .get(reference_sequence_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence ID: {reference_sequence_id}"),
                )
            })?;

        let interval = interval.into();

        let query_bins = reference_sequence
            .query(MIN_SHIFT, DEPTH, interval)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .copied()
            .collect();

        let (start, _) = resolve_interval(interval)?;
        let min_offset = reference_sequence.min_offset(MIN_SHIFT, DEPTH, start);
        let merged_chunks = optimize_chunks(&chunks, min_offset);

        Ok(merged_chunks)
    }
}

impl Default for Index {
    fn default() -> Self {
        Builder::default().build()
    }
}

fn resolve_interval<I>(interval: I) -> io::Result<(Position, Position)>
where
    I: Into<Interval>,
{
    let interval = interval.into();

    let start = interval.start().unwrap_or(Position::MIN);

    if start > MAX_POSITION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid start bound",
        ));
    }

    let end = interval.end().unwrap_or(MAX_POSITION);

    if end > MAX_POSITION {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid end bound",
        ))
    } else {
        Ok((start, end))
    }
}
