//! BAM index bin and fields.

mod builder;

pub(crate) use self::builder::Builder;

use noodles_bgzf::index::Chunk;

use super::DEPTH;

pub(crate) const MAX_ID: u32 = ((1 << (3 * (DEPTH + 1))) - 1) / 7;

pub(crate) const METADATA_ID: u32 = MAX_ID + 1;
pub(crate) const METADATA_CHUNK_COUNT: u32 = 2;

/// A bin in a BAM index reference sequence.
///
/// Bin numbers have an effective range between 0 and 37449, inclusive. An optional pseudo-bin at
/// bin number 37450 holds two pairs of metadata: virtual positions of the start and end of the
/// reference sequence and the number of mapped and unmapped reads in the reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bin {
    id: u32,
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
    pub fn new(id: u32, chunks: Vec<Chunk>) -> Self {
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
    pub fn id(&self) -> u32 {
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
