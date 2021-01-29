//! Tabix index reference sequence and fields.

pub mod bin;
mod builder;
pub mod metadata;

pub use self::{bin::Bin, metadata::Metadata};

pub(crate) use self::builder::Builder;

use std::convert::TryFrom;

use noodles_bgzf as bgzf;

const WINDOW_SIZE: i32 = 16384;

/// A tabix index reference sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReferenceSequence {
    bins: Vec<Bin>,
    intervals: Vec<bgzf::VirtualPosition>,
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
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// ```
    pub fn new(bins: Vec<Bin>, intervals: Vec<bgzf::VirtualPosition>) -> Self {
        Self { bins, intervals }
    }

    /// Returns the list of bins in the reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix::index::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
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
    /// let reference_sequence = ReferenceSequence::new(Vec::new(), Vec::new());
    /// assert!(reference_sequence.intervals().is_empty());
    /// ```
    pub fn intervals(&self) -> &[bgzf::VirtualPosition] {
        &self.intervals
    }

    /// Returns metadata for this reference sequence.
    ///
    /// Metadata is parsed from optional bin 37450.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use noodles_tabix::index::{
    ///     reference_sequence::{bin::Chunk, Bin, Metadata},
    ///     ReferenceSequence,
    /// };
    ///
    /// let reference_sequence = ReferenceSequence::new(
    ///     vec![
    ///         Bin::new(0, Vec::new()),
    ///         Bin::new(37450, vec![
    ///             Chunk::new(bgzf::VirtualPosition::from(610), bgzf::VirtualPosition::from(1597)),
    ///             Chunk::new(bgzf::VirtualPosition::from(55), bgzf::VirtualPosition::from(0)),
    ///         ]),
    ///     ],
    ///     Vec::new(),
    /// );
    ///
    /// assert_eq!(reference_sequence.metadata(), Some(Metadata::new(
    ///     bgzf::VirtualPosition::from(610),
    ///     bgzf::VirtualPosition::from(1597),
    ///     55,
    ///     0,
    /// )));
    /// ```
    pub fn metadata(&self) -> Option<Metadata> {
        self.bins()
            .iter()
            .find(|b| b.id() == metadata::MAGIC_NUMBER)
            .and_then(|bin| Metadata::try_from(bin).ok())
    }
}
