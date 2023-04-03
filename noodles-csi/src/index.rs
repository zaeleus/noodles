//! Coordinate-sorted index and fields.

mod builder;
pub mod header;
mod indexer;
pub mod reference_sequence;

pub use self::{
    builder::Builder, header::Header, indexer::Indexer, reference_sequence::ReferenceSequence,
};

use std::io;

use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Position};

use super::index::reference_sequence::bin::Chunk;

/// A coordinate-sorted index (CSI).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Index {
    min_shift: u8,
    depth: u8,
    header: Option<Header>,
    reference_sequences: Vec<ReferenceSequence>,
    n_no_coor: Option<u64>,
}

impl Index {
    /// Returns a builder to create an index from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let builder = csi::Index::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the number of bits for the minimum interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert_eq!(index.min_shift(), 14);
    /// ```
    pub fn min_shift(&self) -> u8 {
        self.min_shift
    }

    /// Returns the depth of the binning index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert_eq!(index.depth(), 5);
    /// ```
    pub fn depth(&self) -> u8 {
        self.depth
    }

    /// Returns the tabix header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert!(index.header().is_none());
    /// ```
    pub fn header(&self) -> Option<&Header> {
        self.header.as_ref()
    }

    /// Returns a list of indexed reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert!(index.reference_sequences().is_empty());
    /// ```
    pub fn reference_sequences(&self) -> &[ReferenceSequence] {
        &self.reference_sequences
    }

    /// Returns the number of unplaced, unmapped records in the associated file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert!(index.unplaced_unmapped_record_count().is_none());
    /// ```
    pub fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        self.n_no_coor
    }

    /// Returns the chunks that overlap with the given region.
    pub fn query<I>(&self, reference_sequence_id: usize, interval: I) -> io::Result<Vec<Chunk>>
    where
        I: Into<Interval>,
    {
        use super::binning_index::optimize_chunks;

        let interval = interval.into();

        let reference_sequence = self
            .reference_sequences()
            .get(reference_sequence_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence ID: {reference_sequence_id}"),
                )
            })?;

        let query_bins = reference_sequence
            .query(self.min_shift(), self.depth(), interval)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .copied()
            .collect();

        let (start, _) = resolve_interval(self.min_shift(), self.depth(), interval)?;
        let min_offset = reference_sequence.min_offset(self.min_shift(), self.depth(), start);
        let merged_chunks = optimize_chunks(&chunks, min_offset);

        Ok(merged_chunks)
    }

    /// Returns the start position of the first record in the last linear bin.
    ///
    /// This is the closest position to the unplaced, unmapped records, if any, that is available
    /// in an index.
    pub fn first_record_in_last_linear_bin_start_position(&self) -> Option<bgzf::VirtualPosition> {
        self.reference_sequences()
            .iter()
            .rev()
            .find_map(|rs| rs.first_record_in_last_linear_bin_start_position())
    }
}

impl Default for Index {
    fn default() -> Self {
        Self::builder().build()
    }
}

fn resolve_interval<I>(min_shift: u8, depth: u8, interval: I) -> io::Result<(Position, Position)>
where
    I: Into<Interval>,
{
    let interval = interval.into();

    let start = interval.start().unwrap_or(Position::MIN);

    let max_position = ReferenceSequence::max_position(min_shift, depth)?;

    if start > max_position {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid start bound",
        ));
    }

    let end = interval.end().unwrap_or(max_position);

    if end > max_position {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid end bound",
        ))
    } else {
        Ok((start, end))
    }
}
