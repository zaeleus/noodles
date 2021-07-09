//! BAM index reference sequence metadata.

use noodles_bgzf::{index::Chunk, VirtualPosition};

use super::{bin::METADATA_ID, Bin};

/// BAM index reference sequence metadata.
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
    /// use noodles_bgzf as bgzf;
    ///
    /// let metadata = Metadata::new(
    ///     bgzf::VirtualPosition::from(610),
    ///     bgzf::VirtualPosition::from(1597),
    ///     55,
    ///     0,
    /// );
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
    /// use noodles_bgzf as bgzf;
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
    pub fn start_position(&self) -> VirtualPosition {
        self.start_position
    }

    /// Returns the end virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Metadata;
    /// use noodles_bgzf as bgzf;
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
    pub fn end_position(&self) -> VirtualPosition {
        self.end_position
    }

    /// Returns the number of mapped records.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai::index::reference_sequence::Metadata;
    /// use noodles_bgzf as bgzf;;
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
    /// use noodles_bam::bai::index::reference_sequence::Metadata;
    /// use noodles_bgzf as bgzf;;
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

impl From<Metadata> for Bin {
    fn from(metadata: Metadata) -> Self {
        let positions_chunk = Chunk::new(metadata.start_position(), metadata.end_position());

        let counts_chunk = Chunk::new(
            VirtualPosition::from(metadata.mapped_record_count()),
            VirtualPosition::from(metadata.unmapped_record_count()),
        );

        let chunks = vec![positions_chunk, counts_chunk];

        Self::new(METADATA_ID, chunks)
    }
}
