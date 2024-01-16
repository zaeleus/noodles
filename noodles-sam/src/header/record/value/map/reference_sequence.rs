//! SAM header record reference sequence map value.

pub mod alternative_locus;
pub mod alternative_names;
mod builder;
pub mod md5_checksum;
pub mod molecule_topology;
pub mod name;
pub mod tag;

use std::num::NonZeroUsize;

pub(crate) use self::tag::Tag;
pub use self::{
    alternative_locus::AlternativeLocus, alternative_names::AlternativeNames,
    md5_checksum::Md5Checksum, molecule_topology::MoleculeTopology, name::Name,
};

use self::builder::Builder;
use super::{Inner, Map, OtherFields};

/// A SAM header record reference sequence map value.
///
/// The reference sequence describes a sequence a read possibly mapped to. The length is guaranteed
/// to be set.
///
/// A list of reference sequences creates a reference sequence dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    pub(crate) length: NonZeroUsize,
}

impl Inner for ReferenceSequence {
    type StandardTag = tag::Standard;
    type Builder = Builder;
}

impl Map<ReferenceSequence> {
    /// Creates a reference sequence with a length.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn new(length: NonZeroUsize) -> Self {
        Self {
            inner: ReferenceSequence { length },
            other_fields: OtherFields::new(),
        }
    }

    /// Returns the reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert_eq!(usize::from(reference_sequence.length()), 13);
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn length(&self) -> NonZeroUsize {
        self.inner.length
    }

    /// Returns a mutable reference to the reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let length = NonZeroUsize::try_from(13)?;
    /// let mut reference_sequence = Map::<ReferenceSequence>::new(length);
    /// assert_eq!(reference_sequence.length(), length);
    ///
    /// let length = NonZeroUsize::try_from(8)?;
    /// *reference_sequence.length_mut() = length;
    /// assert_eq!(reference_sequence.length(), length);
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn length_mut(&mut self) -> &mut NonZeroUsize {
        &mut self.inner.length
    }
}
