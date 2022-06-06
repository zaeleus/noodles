//! BAM index and fields.

mod builder;
pub mod reference_sequence;

pub use self::{builder::Builder, reference_sequence::ReferenceSequence};

use std::io;

use noodles_core::{region::Interval, Position};
use noodles_csi::{
    binning_index::optimize_chunks, index::reference_sequence::bin::Chunk, BinningIndex,
};

const MIN_SHIFT: u8 = 14;
const DEPTH: u8 = 5;

const MAX_POSITION: Position = match Position::new((1 << (MIN_SHIFT + 3 * DEPTH)) - 1) {
    Some(position) => position,
    None => panic!(),
};

/// A BAM index.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Index {
    reference_sequences: Vec<ReferenceSequence>,
    n_no_coor: Option<u64>,
}

impl Index {
    /// Creates a BAM index builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let builder = bai::Index::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a BAM index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let index = bai::Index::new(Vec::new(), None);
    /// ```
    pub fn new(reference_sequences: Vec<ReferenceSequence>, n_no_coor: Option<u64>) -> Self {
        Self {
            reference_sequences,
            n_no_coor,
        }
    }

    /// Returns the number of unplaced unmapped reads in the associated BAM file.
    ///
    /// An unplaced unmapped read is a read that is has neither a reference sequence ID nor
    /// position set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    ///
    /// let index = bai::Index::default();
    /// assert_eq!(index.unplaced_unmapped_read_count(), None);
    ///
    /// let index = bai::Index::new(Vec::new(), Some(13));
    /// assert_eq!(index.unplaced_unmapped_read_count(), Some(13));
    /// ```
    #[deprecated(
        since = "0.2.0",
        note = "Use `unplaced_unmapped_record_count` instead."
    )]
    pub fn unplaced_unmapped_read_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}

impl BinningIndex for Index {
    type ReferenceSequence = ReferenceSequence;

    /// Returns a list of indexed reference sequences.
    ///
    /// This list is parallel to the reference sequences defined in the associated BAM file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use noodles_csi::BinningIndex;
    /// let index = bai::Index::default();
    /// assert!(index.reference_sequences().is_empty());
    /// ```
    fn reference_sequences(&self) -> &[Self::ReferenceSequence] {
        &self.reference_sequences
    }

    /// Returns the number of unplaced, unmapped records in the associated file.
    ///
    /// An unplaced, unmapped record is one that is has neither a reference sequence ID nor
    /// position set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use noodles_csi::BinningIndex;
    /// let index = bai::Index::default();
    /// assert!(index.unplaced_unmapped_record_count().is_none());
    /// ```
    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        self.n_no_coor
    }

    fn query<I>(&self, reference_sequence_id: usize, interval: I) -> io::Result<Vec<Chunk>>
    where
        I: Into<Interval>,
    {
        let interval = interval.into();

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
            .query(interval)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .copied()
            .collect();

        let (start, _) = resolve_interval(interval)?;
        let min_offset = reference_sequence.min_offset(start);
        let merged_chunks = optimize_chunks(&chunks, min_offset);

        Ok(merged_chunks)
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
