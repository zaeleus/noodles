//! SAM header record reference sequence map value.

pub mod alternative_locus;
pub mod alternative_names;
mod builder;
pub mod md5_checksum;
pub mod molecule_topology;
pub mod name;
mod tag;

use std::{fmt, num::NonZeroUsize};

pub use self::{
    alternative_locus::AlternativeLocus, alternative_names::AlternativeNames,
    md5_checksum::Md5Checksum, molecule_topology::MoleculeTopology, name::Name,
};

use self::{
    builder::Builder,
    tag::{StandardTag, Tag},
};
use super::{Fields, Inner, Map, OtherFields, TryFromFieldsError};

/// A SAM header record reference sequence map value.
///
/// The reference sequence describes a sequence a read possibly mapped to. The length is guaranteed
/// to be set.
///
/// A list of reference sequences creates a reference sequence dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    length: NonZeroUsize,
    alternative_locus: Option<AlternativeLocus>,
    alternative_names: Option<AlternativeNames>,
    assembly_id: Option<String>,
    description: Option<String>,
    md5_checksum: Option<Md5Checksum>,
    species: Option<String>,
    molecule_topology: Option<MoleculeTopology>,
    uri: Option<String>,
}

impl Inner for ReferenceSequence {
    type StandardTag = StandardTag;
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
    pub fn assembly_id(&self) -> Option<&str> {
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
    pub fn description(&self) -> Option<&str> {
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
    pub fn species(&self) -> Option<&str> {
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
    pub fn uri(&self) -> Option<&str> {
        self.inner.uri.as_deref()
    }
}

impl fmt::Display for Map<ReferenceSequence> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "\t{}:{}", tag::LENGTH, self.length())?;

        if let Some(alternative_locus) = self.alternative_locus() {
            write!(f, "\t{}:{alternative_locus}", tag::ALTERNATIVE_LOCUS)?;
        }

        if let Some(alternative_names) = self.alternative_names() {
            write!(f, "\t{}:{alternative_names}", tag::ALTERNATIVE_NAMES)?;
        }

        if let Some(assembly_id) = self.assembly_id() {
            write!(f, "\t{}:{assembly_id}", tag::ASSEMBLY_ID)?;
        }

        if let Some(description) = self.description() {
            write!(f, "\t{}:{description}", tag::DESCRIPTION)?;
        }

        if let Some(md5_checksum) = self.md5_checksum() {
            write!(f, "\t{}:{md5_checksum}", tag::MD5_CHECKSUM)?;
        }

        if let Some(species) = self.species() {
            write!(f, "\t{}:{species}", tag::SPECIES)?;
        }

        if let Some(molecule_topology) = self.molecule_topology() {
            write!(f, "\t{}:{molecule_topology}", tag::MOLECULE_TOPOLOGY)?;
        }

        if let Some(uri) = self.uri() {
            write!(f, "\t{}:{uri}", tag::URI)?;
        }

        super::fmt_display_other_fields(f, self.other_fields())?;

        Ok(())
    }
}

impl TryFrom<Fields> for Map<ReferenceSequence> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut length = None;
        let mut alternative_locus = None;
        let mut alternative_names = None;
        let mut assembly_id = None;
        let mut description = None;
        let mut md5_checksum = None;
        let mut species = None;
        let mut molecule_topology = None;
        let mut uri = None;

        let mut other_fields = OtherFields::new();

        for (key, value) in fields {
            let tag = key.parse().map_err(|_| TryFromFieldsError::InvalidTag)?;

            match tag {
                tag::NAME => return Err(TryFromFieldsError::DuplicateTag),
                tag::LENGTH => {
                    length = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("LN"))?;
                }
                tag::ALTERNATIVE_LOCUS => {
                    alternative_locus = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("AH"))?;
                }
                tag::ALTERNATIVE_NAMES => {
                    alternative_names = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("AN"))?;
                }
                tag::ASSEMBLY_ID => assembly_id = Some(value),
                tag::DESCRIPTION => description = Some(value),
                tag::MD5_CHECKSUM => {
                    md5_checksum = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("M5"))?;
                }
                tag::SPECIES => species = Some(value),
                tag::MOLECULE_TOPOLOGY => {
                    molecule_topology = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("TP"))?;
                }
                tag::URI => uri = Some(value),
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
            }
        }

        let length = length.ok_or(TryFromFieldsError::MissingField("LN"))?;

        Ok(Self {
            inner: ReferenceSequence {
                length,
                alternative_locus,
                alternative_names,
                assembly_id,
                description,
                md5_checksum,
                species,
                molecule_topology,
                uri,
            },
            other_fields,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = Map::<ReferenceSequence>::builder()
            .set_length(NonZeroUsize::try_from(13)?)
            .set_md5_checksum(Md5Checksum::from([
                0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6,
                0x15, 0x34,
            ]))
            .build()?;

        assert_eq!(
            reference_sequence.to_string(),
            "\tLN:13\tM5:d7eba311421bbc9d3ada44709dd61534"
        );

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_duplicate_name() {
        let fields = vec![(String::from("SN"), String::from("sq0"))];

        assert_eq!(
            Map::<ReferenceSequence>::try_from(fields),
            Err(TryFromFieldsError::DuplicateTag)
        );
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_missing_length() {
        assert_eq!(
            Map::<ReferenceSequence>::try_from(vec![]),
            Err(TryFromFieldsError::MissingField("LN"))
        );
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_invalid_length() {
        let fields = vec![(String::from("LN"), String::from("NA"))];

        assert_eq!(
            Map::<ReferenceSequence>::try_from(fields),
            Err(TryFromFieldsError::InvalidValue("LN"))
        );

        let fields = vec![(String::from("LN"), String::from("0"))];

        assert_eq!(
            Map::<ReferenceSequence>::try_from(fields),
            Err(TryFromFieldsError::InvalidValue("LN"))
        );
    }
}
