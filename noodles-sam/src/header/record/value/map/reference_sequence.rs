//! SAM header record reference sequence map value.

pub mod alternative_locus;
pub mod alternative_names;
mod builder;
pub mod md5_checksum;
pub mod molecule_topology;
pub mod name;
mod tag;

use std::{error, fmt, num::NonZeroUsize};

pub use self::{
    alternative_locus::AlternativeLocus, alternative_names::AlternativeNames,
    md5_checksum::Md5Checksum, molecule_topology::MoleculeTopology, name::Name,
};

use self::builder::Builder;
use super::{Fields, Inner, Map, OtherFields, TryFromFieldsError};

type StandardTag = tag::Standard;
type Tag = super::tag::Tag<StandardTag>;

/// A SAM header record reference sequence map value.
///
/// The reference sequence describes a sequence a read possibly mapped to. Both the reference
/// sequence name and length are guaranteed to be set.
///
/// A list of reference sequences creates a reference sequence dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    name: Name,
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

/// An error returned when a SAM header reference sequence fails to construct.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum NewError {
    /// The name is invalid.
    InvalidName,
    /// The length is invalid.
    InvalidLength(usize),
}

impl error::Error for NewError {}

impl fmt::Display for NewError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidName => f.write_str("invalid name"),
            Self::InvalidLength(length) => write!(f, "invalid length: {}", length),
        }
    }
}

impl Map<ReferenceSequence> {
    /// Creates a reference sequence with a name and length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(name: Name, length: usize) -> Result<Self, NewError> {
        let length = NonZeroUsize::new(length).ok_or(NewError::InvalidLength(length))?;

        Ok(Self {
            inner: ReferenceSequence {
                name,
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
        })
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert_eq!(**reference_sequence.name(), "sq0");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn name(&self) -> &Name {
        &self.inner.name
    }

    /// Returns a mutable reference to the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert_eq!(**reference_sequence.name(), "sq0");
    ///
    /// *reference_sequence.name_mut() = "sq1".parse()?;
    /// assert_eq!(**reference_sequence.name(), "sq1");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn name_mut(&mut self) -> &mut Name {
        &mut self.inner.name
    }

    /// Returns the reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert_eq!(usize::from(reference_sequence.length()), 13);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
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
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert_eq!(usize::from(reference_sequence.length()), 13);
    ///
    /// *reference_sequence.length_mut() = NonZeroUsize::try_from(8)?;
    /// assert_eq!(usize::from(reference_sequence.length()), 8);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn length_mut(&mut self) -> &mut NonZeroUsize {
        &mut self.inner.length
    }

    /// Returns the alternative locus.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.alternative_locus().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternative_locus(&self) -> Option<&AlternativeLocus> {
        self.inner.alternative_locus.as_ref()
    }

    /// Returns the alternative names (aliases) of the reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.alternative_names().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternative_names(&self) -> Option<&AlternativeNames> {
        self.inner.alternative_names.as_ref()
    }

    /// Returns the genome assembly ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.assembly_id().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn assembly_id(&self) -> Option<&str> {
        self.inner.assembly_id.as_deref()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.description().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.inner.description.as_deref()
    }

    /// Returns the MD5 checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.md5_checksum().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn md5_checksum(&self) -> Option<Md5Checksum> {
        self.inner.md5_checksum
    }

    /// Returns a mutable reference to the MD5 checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::Md5Checksum, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
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
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.species().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn species(&self) -> Option<&str> {
        self.inner.species.as_deref()
    }

    /// Returns the molecule topology.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.molecule_topology().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn molecule_topology(&self) -> Option<MoleculeTopology> {
        self.inner.molecule_topology
    }

    /// Returns the URI.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    /// let mut reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.uri().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn uri(&self) -> Option<&str> {
        self.inner.uri.as_deref()
    }
}

impl fmt::Display for Map<ReferenceSequence> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SN:{}", self.name())?;
        write!(f, "\tLN:{}", self.length())?;

        if let Some(alternative_locus) = self.alternative_locus() {
            write!(f, "\tAH:{}", alternative_locus)?;
        }

        if let Some(alternative_names) = self.alternative_names() {
            write!(f, "\tAN:{}", alternative_names)?;
        }

        if let Some(assembly_id) = self.assembly_id() {
            write!(f, "\tAS:{}", assembly_id)?;
        }

        if let Some(description) = self.description() {
            write!(f, "\tDS:{}", description)?;
        }

        if let Some(md5_checksum) = self.md5_checksum() {
            write!(f, "\tM5:{}", md5_checksum)?;
        }

        if let Some(species) = self.species() {
            write!(f, "\tSP:{}", species)?;
        }

        if let Some(molecule_topology) = self.molecule_topology() {
            write!(f, "\tTP:{}", molecule_topology)?;
        }

        if let Some(uri) = self.uri() {
            write!(f, "\tUR:{}", uri)?;
        }

        super::fmt_display_other_fields(f, self.other_fields())?;

        Ok(())
    }
}

impl TryFrom<Fields> for Map<ReferenceSequence> {
    type Error = TryFromFieldsError;

    fn try_from(fields: Fields) -> Result<Self, Self::Error> {
        let mut other_fields = super::init_other_fields(fields.len());

        let mut name = None;
        let mut length = None;
        let mut alternative_locus = None;
        let mut alternative_names = None;
        let mut assembly_id = None;
        let mut description = None;
        let mut md5_checksum = None;
        let mut species = None;
        let mut molecule_topology = None;
        let mut uri = None;

        for (key, value) in fields {
            let tag = key.parse().map_err(|_| TryFromFieldsError::InvalidTag)?;

            match tag {
                Tag::Standard(StandardTag::Name) => {
                    name = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("SN"))?;
                }
                Tag::Standard(StandardTag::Length) => {
                    length = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("LN"))?;
                }
                Tag::Standard(StandardTag::AlternativeLocus) => {
                    alternative_locus = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("AH"))?;
                }
                Tag::Standard(StandardTag::AlternativeNames) => {
                    alternative_names = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("AN"))?;
                }
                Tag::Standard(StandardTag::AssemblyId) => assembly_id = Some(value),
                Tag::Standard(StandardTag::Description) => description = Some(value),
                Tag::Standard(StandardTag::Md5Checksum) => {
                    md5_checksum = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("M5"))?;
                }
                Tag::Standard(StandardTag::Species) => species = Some(value),
                Tag::Standard(StandardTag::MoleculeTopology) => {
                    molecule_topology = value
                        .parse()
                        .map(Some)
                        .map_err(|_| TryFromFieldsError::InvalidValue("TP"))?;
                }
                Tag::Standard(StandardTag::Uri) => uri = Some(value),
                Tag::Other(t) => super::insert_other_field(&mut other_fields, t, value)?,
            }
        }

        let name = name.ok_or(TryFromFieldsError::MissingField("SN"))?;
        let length = length.ok_or(TryFromFieldsError::MissingField("LN"))?;

        Ok(Self {
            inner: ReferenceSequence {
                name,
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
    fn test_new() -> Result<(), name::ParseError> {
        assert_eq!(
            Map::<ReferenceSequence>::new("sq0".parse()?, 0),
            Err(NewError::InvalidLength(0))
        );

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = Map::<ReferenceSequence>::builder()
            .set_name("sq0".parse()?)
            .set_length(13)
            .set_md5_checksum(Md5Checksum::from([
                0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6,
                0x15, 0x34,
            ]))
            .build()?;

        assert_eq!(
            reference_sequence.to_string(),
            "SN:sq0\tLN:13\tM5:d7eba311421bbc9d3ada44709dd61534"
        );

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_missing_name(
    ) -> Result<(), TryFromFieldsError> {
        let fields = vec![(String::from("LN"), String::from("8"))];

        assert_eq!(
            Map::<ReferenceSequence>::try_from(fields),
            Err(TryFromFieldsError::MissingField("SN"))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_missing_length(
    ) -> Result<(), TryFromFieldsError> {
        let fields = vec![(String::from("SN"), String::from("sq0"))];

        assert_eq!(
            Map::<ReferenceSequence>::try_from(fields),
            Err(TryFromFieldsError::MissingField("LN"))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_invalid_name(
    ) -> Result<(), TryFromFieldsError> {
        let fields = vec![(String::from("SN"), String::from("*"))];

        assert_eq!(
            Map::<ReferenceSequence>::try_from(fields),
            Err(TryFromFieldsError::InvalidValue("SN"))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_invalid_length(
    ) -> Result<(), TryFromFieldsError> {
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

        Ok(())
    }
}
