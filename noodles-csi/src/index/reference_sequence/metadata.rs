//! CSI reference sequence metadata.

use std::{convert::TryFrom, error, fmt};

use noodles_bgzf as bgzf;

use super::Bin;

/// CSI reference sequence metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Metadata {
    start_position: bgzf::VirtualPosition,
    end_position: bgzf::VirtualPosition,
    mapped_record_count: u64,
    unmapped_record_count: u64,
}

impl Metadata {
    /// Creates CSI reference sequence metadata.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Metadata;
    ///
    /// let metadata = Metadata::new(
    ///     bgzf::VirtualPosition::from(610),
    ///     bgzf::VirtualPosition::from(1597),
    ///     55,
    ///     0,
    /// );
    /// ```
    pub fn new(
        start_position: bgzf::VirtualPosition,
        end_position: bgzf::VirtualPosition,
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
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Metadata;
    ///
    /// let metadata = Metadata::new(
    ///     bgzf::VirtualPosition::from(610),
    ///     bgzf::VirtualPosition::from(1597),
    ///     55,
    ///     0,
    /// );
    ///
    /// assert_eq!(metadata.start_position(), bgzf::VirtualPosition::from(610));
    /// ```
    pub fn start_position(&self) -> bgzf::VirtualPosition {
        self.start_position
    }

    /// Returns the end virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Metadata;
    ///
    /// let metadata = Metadata::new(
    ///     bgzf::VirtualPosition::from(610),
    ///     bgzf::VirtualPosition::from(1597),
    ///     55,
    ///     0,
    /// );
    ///
    /// assert_eq!(metadata.end_position(), bgzf::VirtualPosition::from(1597));
    /// ```
    pub fn end_position(&self) -> bgzf::VirtualPosition {
        self.end_position
    }

    /// Returns the number of mapped records.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Metadata;
    ///
    /// let metadata = Metadata::new(
    ///     bgzf::VirtualPosition::from(610),
    ///     bgzf::VirtualPosition::from(1597),
    ///     55,
    ///     0,
    /// );
    ///
    /// assert_eq!(metadata.mapped_record_count(), 55);
    /// ```
    pub fn mapped_record_count(&self) -> u64 {
        self.mapped_record_count
    }

    /// Returns the number of unmapped records.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::Metadata;
    ///
    /// let metadata = Metadata::new(
    ///     bgzf::VirtualPosition::from(610),
    ///     bgzf::VirtualPosition::from(1597),
    ///     55,
    ///     0,
    /// );
    ///
    /// assert_eq!(metadata.unmapped_record_count(), 0);
    /// ```
    pub fn unmapped_record_count(&self) -> u64 {
        self.unmapped_record_count
    }
}

/// An error returned when a raw bin fails to convert to metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromBinError {
    /// The positions chunk is missing.
    MissingPositionsChunk,
    /// The counts chunk is missing.
    MissingCountsChunk,
}

impl error::Error for TryFromBinError {}

impl fmt::Display for TryFromBinError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPositionsChunk => f.write_str("missing positions chunk"),
            Self::MissingCountsChunk => f.write_str("missing counts chunk"),
        }
    }
}

impl TryFrom<&Bin> for Metadata {
    type Error = TryFromBinError;

    fn try_from(bin: &Bin) -> Result<Self, Self::Error> {
        let mut chunks_iter = bin.chunks().iter();

        let (ref_beg, ref_end) = chunks_iter
            .next()
            .map(|c| (c.start(), c.end()))
            .ok_or(TryFromBinError::MissingPositionsChunk)?;

        let (n_mapped, n_unmapped) = chunks_iter
            .next()
            .map(|c| (u64::from(c.start()), u64::from(c.end())))
            .ok_or(TryFromBinError::MissingCountsChunk)?;

        Ok(Self::new(ref_beg, ref_end, n_mapped, n_unmapped))
    }
}
