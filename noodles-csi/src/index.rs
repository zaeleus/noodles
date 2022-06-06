//! Coordinate-sorted index and fields.

mod builder;
pub mod reference_sequence;

pub use self::{builder::Builder, reference_sequence::ReferenceSequence};

use std::io;

use noodles_core::{region::Interval, Position};

use super::{index::reference_sequence::bin::Chunk, BinningIndex};

/// A coordinate-sorted index (CSI).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Index {
    min_shift: u8,
    depth: u8,
    aux: Vec<u8>,
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

    /// Returns the auxiliary data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert!(index.aux().is_empty());
    /// ```
    pub fn aux(&self) -> &[u8] {
        &self.aux
    }

    /// Returns the number of unmapped records in the associated file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// assert!(index.unmapped_read_count().is_none());
    /// ```
    #[deprecated(
        since = "0.2.0",
        note = "Use `unplaced_unmapped_record_count` instead."
    )]
    pub fn unmapped_read_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}

impl BinningIndex for Index {
    type ReferenceSequence = ReferenceSequence;

    /// Returns a list of indexed reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi::{self as csi, BinningIndex};
    /// let index = csi::Index::default();
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
    /// use noodles_csi::{self as csi, BinningIndex};
    /// let index = csi::Index::default();
    /// assert!(index.unplaced_unmapped_record_count().is_none());
    /// ```
    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        self.n_no_coor
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
                    format!("invalid reference sequence ID: {}", reference_sequence_id),
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

        Ok(chunks)
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
