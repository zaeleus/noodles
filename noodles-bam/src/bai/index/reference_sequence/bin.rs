//! BAM index bin and fields.

mod builder;

pub(crate) use self::builder::Builder;

use noodles_csi::index::reference_sequence::bin::Chunk;

use crate::bai::index::DEPTH;

pub(crate) const MAX_ID: usize = ((1 << (3 * (DEPTH + 1))) - 1) / 7;

pub(crate) const METADATA_ID: usize = MAX_ID + 1;
pub(crate) const METADATA_CHUNK_COUNT: usize = 2;

/// A bin in a BAM index reference sequence.
///
/// Bin numbers have an effective range between 0 and 37449, inclusive. An optional pseudo-bin at
/// bin number 37450 holds two pairs of metadata: virtual positions of the start and end of the
/// reference sequence and the number of mapped and unmapped reads in the reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bin {
    id: usize,
    chunks: Vec<Chunk>,
}

impl Bin {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a BAM index reference sequence bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// ```
    pub fn new(id: usize, chunks: Vec<Chunk>) -> Self {
        Self { id, chunks }
    }

    /// Returns the bin ID.
    ///
    /// This is also called the bin number.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// assert_eq!(bin.id(), 10946);
    /// ```
    pub fn id(&self) -> usize {
        self.id
    }

    /// Returns the list of chunks in the bin.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Bin;
    /// let bin = Bin::new(10946, Vec::new());
    /// assert!(bin.chunks().is_empty());
    /// ```
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }
}
