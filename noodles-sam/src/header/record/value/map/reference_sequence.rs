//! SAM header record reference sequence map value.

pub mod alternative_locus;
pub mod alternative_names;
mod builder;
pub mod md5_checksum;
pub mod molecule_topology;
pub mod name;
pub(crate) mod tag;

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
    pub(crate) alternative_locus: Option<AlternativeLocus>,
    pub(crate) alternative_names: Option<AlternativeNames>,
    pub(crate) assembly_id: Option<Vec<u8>>,
    pub(crate) description: Option<Vec<u8>>,
    pub(crate) md5_checksum: Option<Md5Checksum>,
    pub(crate) species: Option<Vec<u8>>,
    pub(crate) molecule_topology: Option<MoleculeTopology>,
    pub(crate) uri: Option<Vec<u8>>,
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
            inner: ReferenceSequence {
                length,
                alternative_locus: None,
                alternative_names: None,
                assembly_id: None,
                description: None,
                md5_checksum: None,
                species: None,
                molecule_topology: None,
                uri: None,
            },
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

    /// Returns the alternative locus.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.alternative_locus().is_none());
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn alternative_locus(&self) -> Option<&AlternativeLocus> {
        self.inner.alternative_locus.as_ref()
    }

    /// Returns the alternative names (aliases) of the reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.alternative_names().is_none());
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn alternative_names(&self) -> Option<&AlternativeNames> {
        self.inner.alternative_names.as_ref()
    }

    /// Returns the genome assembly ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.assembly_id().is_none());
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn assembly_id(&self) -> Option<&[u8]> {
        self.inner.assembly_id.as_deref()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.description().is_none());
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn description(&self) -> Option<&[u8]> {
        self.inner.description.as_deref()
    }

    /// Returns the MD5 checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.md5_checksum().is_none());
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn md5_checksum(&self) -> Option<Md5Checksum> {
        self.inner.md5_checksum
    }

    /// Returns a mutable reference to the MD5 checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::Md5Checksum, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let mut reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.md5_checksum().is_none());
    ///
    /// let checksum: Md5Checksum = "d7eba311421bbc9d3ada44709dd61534".parse()?;
    /// *reference_sequence.md5_checksum_mut() = Some(checksum);
    /// assert_eq!(reference_sequence.md5_checksum(), Some(checksum));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn md5_checksum_mut(&mut self) -> &mut Option<Md5Checksum> {
        &mut self.inner.md5_checksum
    }

    /// Returns the species.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.species().is_none());
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn species(&self) -> Option<&[u8]> {
        self.inner.species.as_deref()
    }

    /// Returns the molecule topology.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.molecule_topology().is_none());
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn molecule_topology(&self) -> Option<MoleculeTopology> {
        self.inner.molecule_topology
    }

    /// Returns the URI.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?);
    /// assert!(reference_sequence.uri().is_none());
    /// # Ok::<_, std::num::TryFromIntError>(())
    /// ```
    pub fn uri(&self) -> Option<&[u8]> {
        self.inner.uri.as_deref()
    }
}
