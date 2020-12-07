//! BAM index reference sequence metadata.

use std::{convert::TryFrom, error, fmt};

use noodles_bgzf::VirtualPosition;

use super::{bin::Chunk, Bin};

pub(crate) const MAGIC_NUMBER: u32 = 37450;

/// Metadatadata in a BAM index reference sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Metadata {
    start_position: VirtualPosition,
    end_position: VirtualPosition,
    mapped_record_count: u64,
    unmapped_record_count: u64,
}

impl Metadata {
    /// Creates reference sequence metadata.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Metadata;
    /// use noodles_bgzf::VirtualPosition;
    /// let meta = Metadata::new(VirtualPosition::from(610), VirtualPosition::from(1597), 55, 0);
    /// ```
    pub fn new(
        start_position: VirtualPosition,
        end_position: VirtualPosition,
        mapped_record_count: u64,
        unmapped_record_count: u64,
    ) -> Self {
        Self {
            start_position,
            end_position,
            mapped_record_count,
            unmapped_record_count,
        }
    }

    /// Returns the start virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Metadata;
    /// use noodles_bgzf::VirtualPosition;
    /// let meta = Metadata::new(VirtualPosition::from(610), VirtualPosition::from(1597), 55, 0);
    /// assert_eq!(meta.start_position(), VirtualPosition::from(610));
    /// ```
    pub fn start_position(&self) -> VirtualPosition {
        self.start_position
    }

    /// Returns the end virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Metadata;
    /// use noodles_bgzf::VirtualPosition;
    /// let meta = Metadata::new(VirtualPosition::from(610), VirtualPosition::from(1597), 55, 0);
    /// assert_eq!(meta.end_position(), VirtualPosition::from(1597));
    /// ```
    pub fn end_position(&self) -> VirtualPosition {
        self.end_position
    }

    /// Returns the number of mapped records.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Metadata;
    /// use noodles_bgzf::VirtualPosition;
    /// let meta = Metadata::new(VirtualPosition::from(610), VirtualPosition::from(1597), 55, 0);
    /// assert_eq!(meta.mapped_record_count(), 55);
    /// ```
    pub fn mapped_record_count(&self) -> u64 {
        self.mapped_record_count
    }

    /// Returns the number of unmapped records.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Metadata;
    /// use noodles_bgzf::VirtualPosition;
    /// let meta = Metadata::new(VirtualPosition::from(610), VirtualPosition::from(1597), 55, 0);
    /// assert_eq!(meta.unmapped_record_count(), 0);
    /// ```
    pub fn unmapped_record_count(&self) -> u64 {
        self.unmapped_record_count
    }
}

/// An error returned when a raw bin fails to convert to metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromBinError {
    /// The bin number is invalid.
    InvalidMagicNumber(u32),
    /// The positions chunk is missing.
    MissingPositionsChunk,
    /// The counts chunk is missing.
    MissingCountsChunk,
}

impl error::Error for TryFromBinError {}

impl fmt::Display for TryFromBinError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMagicNumber(n) => write!(
                f,
                "invalid magic number: expected {}, got {}",
                MAGIC_NUMBER, n
            ),
            Self::MissingPositionsChunk => f.write_str("missing positions chunk"),
            Self::MissingCountsChunk => f.write_str("missing counts chunk"),
        }
    }
}

impl TryFrom<&Bin> for Metadata {
    type Error = TryFromBinError;

    fn try_from(bin: &Bin) -> Result<Self, Self::Error> {
        if bin.id() != MAGIC_NUMBER {
            return Err(TryFromBinError::InvalidMagicNumber(bin.id()));
        }

        let mut chunks_iter = bin.chunks().iter();

        let (ref_beg, ref_end) = chunks_iter
            .next()
            .map(|c| (c.start(), c.end()))
            .ok_or(TryFromBinError::MissingPositionsChunk)?;

        let (n_mapped, n_unmapped) = chunks_iter
            .next()
            .map(|c| (u64::from(c.start()), u64::from(c.end())))
            .ok_or(TryFromBinError::MissingCountsChunk)?;

        Ok(Self {
            start_position: ref_beg,
            end_position: ref_end,
            mapped_record_count: n_mapped,
            unmapped_record_count: n_unmapped,
        })
    }
}

impl From<Metadata> for Bin {
    fn from(metadata: Metadata) -> Self {
        let positions_chunk = Chunk::new(metadata.start_position(), metadata.end_position());

        let counts_chunk = Chunk::new(
            VirtualPosition::from(metadata.mapped_record_count()),
            VirtualPosition::from(metadata.unmapped_record_count()),
        );

        let chunks = vec![positions_chunk, counts_chunk];

        Self::new(MAGIC_NUMBER, chunks)
    }
}
