//! Tabix index reference sequence and fields.

pub mod bin;
mod builder;

pub use self::bin::Bin;

pub(crate) use self::builder::Builder;

use noodles_bgzf as bgzf;

const WINDOW_SIZE: u32 = 16384;

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
}
