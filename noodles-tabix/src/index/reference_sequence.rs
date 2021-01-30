//! Tabix index reference sequence and fields.

pub mod bin;
mod builder;
pub mod metadata;

pub use self::{bin::Bin, metadata::Metadata};

pub(crate) use self::builder::Builder;

use noodles_bgzf as bgzf;

const WINDOW_SIZE: i32 = 16384;

/// A tabix index reference sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReferenceSequence {
    bins: Vec<Bin>,
    intervals: Vec<bgzf::VirtualPosition>,
    metadata: Option<Metadata>,
}

impl ReferenceSequence {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a tabix index reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// ```
    pub fn new(
        bins: Vec<Bin>,
        intervals: Vec<bgzf::VirtualPosition>,
        metadata: Option<Metadata>,
    ) -> Self {
        Self {
            bins,
            intervals,
            metadata,
        }
    }

    /// Returns the list of bins in the reference sequence.
    ///
    /// This list does include the metadata pseudo-bin (bin 37450). Use [`Self::metadata`] instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert!(reference_sequence.bins().is_empty());
    /// ```
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }

    /// Returns the list of 16 kbp intervals that make up the linear index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert!(reference_sequence.intervals().is_empty());
    /// ```
    pub fn intervals(&self) -> &[bgzf::VirtualPosition] {
        &self.intervals
    }

    /// Returns metadata for this reference sequence.
    ///
    /// Metadata is parsed from the optional pseudo-bin 37450.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::VirtualPosition;
    /// use noodles_tabix::index::{reference_sequence::Metadata, ReferenceSequence};
    ///
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new(), None);
    /// assert!(reference_sequence.metadata().is_none());
    ///
    /// let reference_sequence = ReferenceSequence::new(
    ///     Vec::new(),
    ///     Vec::new(),
    ///     Some(Metadata::new(VirtualPosition::from(610), VirtualPosition::from(1597), 55, 0))
    /// );
    /// assert!(reference_sequence.metadata().is_some());
    /// ```
    pub fn metadata(&self) -> Option<&Metadata> {
        self.metadata.as_ref()
    }
}
