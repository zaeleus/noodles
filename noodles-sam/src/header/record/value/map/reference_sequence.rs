//! SAM header record reference sequence map value.

pub mod alternative_locus;
pub mod alternative_names;
mod builder;
pub mod md5_checksum;
pub mod molecule_topology;
pub mod name;
mod tag;

use std::{
    error, fmt,
    num::{self, NonZeroUsize},
};

pub use self::{
    alternative_locus::AlternativeLocus, alternative_names::AlternativeNames,
    md5_checksum::Md5Checksum, molecule_topology::MoleculeTopology, name::Name,
};

use self::{
    builder::Builder,
    tag::{StandardTag, Tag},
};
use super::{Fields, Inner, Map, OtherFields};

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

/// An error returned when a raw header reference sequence record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A field is missing.
    MissingField(Tag),
    /// A tag is invalid.
    InvalidTag(super::tag::ParseError),
    /// A tag is duplicated.
    DuplicateTag(Tag),
    /// The length is invalid.
    InvalidLength(num::ParseIntError),
    /// The alternative locus is invalid.
    InvalidAlternativeLocus(alternative_locus::ParseError),
    /// The alternative names is invalid.
    InvalidAlternativeNames(alternative_names::ParseError),
    /// The MD5 checksum is invalid.
    InvalidMd5Checksum(md5_checksum::ParseError),
    /// The molecule topology is invalid.
    InvalidMoleculeToplogy(molecule_topology::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidTag(e) => Some(e),
            Self::InvalidLength(e) => Some(e),
            Self::InvalidAlternativeLocus(e) => Some(e),
            Self::InvalidAlternativeNames(e) => Some(e),
            Self::InvalidMd5Checksum(e) => Some(e),
            Self::InvalidMoleculeToplogy(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField(tag) => write!(f, "missing field: {tag}"),
            Self::InvalidTag(_) => write!(f, "invalid tag"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
            Self::InvalidAlternativeLocus(_) => write!(f, "invalid alternative locus"),
            Self::InvalidAlternativeNames(_) => write!(f, "invalid alternative names"),
            Self::InvalidMd5Checksum(_) => write!(f, "invalid MD5 checksum"),
            Self::InvalidMoleculeToplogy(_) => write!(f, "invalid molecule topology"),
        }
    }
}

impl TryFrom<Fields> for Map<ReferenceSequence> {
    type Error = ParseError;

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
            let tag = key.parse().map_err(ParseError::InvalidTag)?;

            match tag {
                tag::NAME => return Err(ParseError::DuplicateTag(tag::NAME)),
                tag::LENGTH => length = parse_length(&value).map(Some)?,
                tag::ALTERNATIVE_LOCUS => {
                    alternative_locus = parse_alternative_locus(&value).map(Some)?
                }
                tag::ALTERNATIVE_NAMES => {
                    alternative_names = parse_alternative_names(&value).map(Some)?
                }
                tag::ASSEMBLY_ID => assembly_id = Some(value),
                tag::DESCRIPTION => description = Some(value),
                tag::MD5_CHECKSUM => md5_checksum = parse_md5_checksum(&value).map(Some)?,
                tag::SPECIES => species = Some(value),
                tag::MOLECULE_TOPOLOGY => {
                    molecule_topology = parse_molecule_topology(&value).map(Some)?
                }
                tag::URI => uri = Some(value),
                Tag::Other(t) => try_insert(&mut other_fields, t, value)?,
            }
        }

        let length = length.ok_or(ParseError::MissingField(tag::LENGTH))?;

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

fn parse_length(s: &str) -> Result<NonZeroUsize, ParseError> {
    s.parse().map_err(ParseError::InvalidLength)
}

fn parse_alternative_locus(s: &str) -> Result<AlternativeLocus, ParseError> {
    s.parse().map_err(ParseError::InvalidAlternativeLocus)
}

fn parse_alternative_names(s: &str) -> Result<AlternativeNames, ParseError> {
    s.parse().map_err(ParseError::InvalidAlternativeNames)
}

fn parse_md5_checksum(s: &str) -> Result<Md5Checksum, ParseError> {
    s.parse().map_err(ParseError::InvalidMd5Checksum)
}

fn parse_molecule_topology(s: &str) -> Result<MoleculeTopology, ParseError> {
    s.parse().map_err(ParseError::InvalidMoleculeToplogy)
}

fn try_insert(
    other_fields: &mut OtherFields<StandardTag>,
    tag: super::tag::Other<StandardTag>,
    value: String,
) -> Result<(), ParseError> {
    use indexmap::map::Entry;

    match other_fields.entry(tag) {
        Entry::Vacant(entry) => {
            entry.insert(value);
            Ok(())
        }
        Entry::Occupied(entry) => {
            let (t, _) = entry.remove_entry();
            Err(ParseError::DuplicateTag(Tag::Other(t)))
        }
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
            Err(ParseError::DuplicateTag(tag::NAME))
        );
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_missing_length() {
        assert_eq!(
            Map::<ReferenceSequence>::try_from(vec![]),
            Err(ParseError::MissingField(tag::LENGTH))
        );
    }

    #[test]
    fn test_try_from_fields_for_map_reference_sequence_with_invalid_length() {
        let fields = vec![(String::from("LN"), String::from("NA"))];

        assert!(matches!(
            Map::<ReferenceSequence>::try_from(fields),
            Err(ParseError::InvalidLength(_))
        ));

        let fields = vec![(String::from("LN"), String::from("0"))];

        assert!(matches!(
            Map::<ReferenceSequence>::try_from(fields),
            Err(ParseError::InvalidLength(_))
        ));
    }
}
