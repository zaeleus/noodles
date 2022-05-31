//! SAM header reference sequence and fields.

pub mod alternative_locus;
pub mod alternative_names;
pub mod builder;
pub mod md5_checksum;
pub mod molecule_topology;
pub mod name;
pub mod tag;

use std::{collections::HashMap, error, fmt, num::NonZeroUsize};

pub use self::{
    alternative_locus::AlternativeLocus, alternative_names::AlternativeNames, builder::Builder,
    md5_checksum::Md5Checksum, molecule_topology::MoleculeTopology, name::Name, tag::Tag,
};

use super::{
    record::{self, value::Fields},
    Record,
};

/// A SAM header reference sequence.
///
/// The reference sequence describes a sequence a read possibly mapped to. Both the reference
/// sequence name and length are guaranteed to be set.
///
/// A list of reference sequences creates a reference sequence dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    name: Name,
    len: NonZeroUsize,
    alternative_locus: Option<AlternativeLocus>,
    alternative_names: Option<AlternativeNames>,
    assembly_id: Option<String>,
    description: Option<String>,
    md5_checksum: Option<Md5Checksum>,
    species: Option<String>,
    molecule_topology: Option<MoleculeTopology>,
    uri: Option<String>,
    fields: HashMap<Tag, String>,
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
            Self::InvalidLength(len) => write!(f, "invalid length: {}", len),
        }
    }
}

#[allow(clippy::len_without_is_empty)]
impl ReferenceSequence {
    /// Creates a SAM header reference sequence builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let builder = ReferenceSequence::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a reference sequence with a name and length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(name: Name, len: usize) -> Result<Self, NewError> {
        let len = NonZeroUsize::new(len).ok_or(NewError::InvalidLength(len))?;

        Ok(Self {
            name,
            len,
            alternative_locus: None,
            alternative_names: None,
            assembly_id: None,
            description: None,
            md5_checksum: None,
            species: None,
            molecule_topology: None,
            uri: None,
            fields: HashMap::new(),
        })
    }

    /// Returns the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert_eq!(**reference_sequence.name(), "sq0");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn name(&self) -> &Name {
        &self.name
    }

    /// Returns a mutable reference to the reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::Name, ReferenceSequence};
    ///
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert_eq!(**reference_sequence.name(), "sq0");
    ///
    /// *reference_sequence.name_mut() = "sq1".parse()?;
    /// assert_eq!(**reference_sequence.name(), "sq1");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn name_mut(&mut self) -> &mut Name {
        &mut self.name
    }

    /// Returns the reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert_eq!(usize::from(reference_sequence.len()), 13);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn len(&self) -> NonZeroUsize {
        self.len
    }

    /// Returns a mutable reference to the reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert_eq!(usize::from(reference_sequence.len()), 13);
    ///
    /// *reference_sequence.len_mut() = NonZeroUsize::try_from(8)?;
    /// assert_eq!(usize::from(reference_sequence.len()), 8);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn len_mut(&mut self) -> &mut NonZeroUsize {
        &mut self.len
    }

    /// Returns the alternative locus.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.alternative_locus().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternative_locus(&self) -> Option<&AlternativeLocus> {
        self.alternative_locus.as_ref()
    }

    /// Returns the alternative names (aliases) of the reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.alternative_names().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn alternative_names(&self) -> Option<&AlternativeNames> {
        self.alternative_names.as_ref()
    }

    /// Returns the genome assembly ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.assembly_id().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn assembly_id(&self) -> Option<&str> {
        self.assembly_id.as_deref()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.description().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Returns the MD5 checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.md5_checksum().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn md5_checksum(&self) -> Option<Md5Checksum> {
        self.md5_checksum
    }

    /// Returns a mutable reference to the MD5 checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::Md5Checksum, ReferenceSequence};
    ///
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.md5_checksum().is_none());
    ///
    /// let checksum: Md5Checksum = "d7eba311421bbc9d3ada44709dd61534".parse()?;
    /// *reference_sequence.md5_checksum_mut() = Some(checksum);
    /// assert_eq!(reference_sequence.md5_checksum(), Some(checksum));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn md5_checksum_mut(&mut self) -> &mut Option<Md5Checksum> {
        &mut self.md5_checksum
    }

    /// Returns the species.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.species().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn species(&self) -> Option<&str> {
        self.species.as_deref()
    }

    /// Returns the molecule topology.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.molecule_topology().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn molecule_topology(&self) -> Option<MoleculeTopology> {
        self.molecule_topology
    }

    /// Returns the URI.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    /// let mut reference_sequence = ReferenceSequence::new("sq0".parse()?, 13)?;
    /// assert!(reference_sequence.uri().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn uri(&self) -> Option<&str> {
        self.uri.as_deref()
    }

    /// Returns the raw fields of the reference sequence.
    ///
    /// This includes any field that is not specially handled by the structure itself. For example,
    /// this will not include the name and length fields, as they are parsed and available as
    /// [`Self::name`] and [`Self::len`], respectively.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::{self, Tag}, ReferenceSequence};
    ///
    /// let name: reference_sequence::Name = "sq0".parse()?;
    /// let zn = Tag::Other([b'z', b'n']);
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name(name.clone())
    ///     .set_length(13)
    ///     .insert(zn, String::from("noodles"))
    ///     .build()?;
    ///
    /// let fields = reference_sequence.fields();
    /// assert_eq!(fields.len(), 1);
    /// assert_eq!(fields.get(&zn), Some(&String::from("noodles")));
    ///
    /// assert_eq!(fields.get(&Tag::Name), None);
    /// assert_eq!(reference_sequence.name(), &name);
    ///
    /// assert_eq!(fields.get(&Tag::Length), None);
    /// assert_eq!(usize::from(reference_sequence.len()), 13);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn fields(&self) -> &HashMap<Tag, String> {
        &self.fields
    }
}

impl fmt::Display for ReferenceSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::ReferenceSequence)?;
        write!(f, "\t{}:{}", Tag::Name, self.name())?;
        write!(f, "\t{}:{}", Tag::Length, self.len())?;

        if let Some(alternative_locus) = self.alternative_locus() {
            write!(f, "\t{}:{}", Tag::AlternativeLocus, alternative_locus)?;
        }

        if let Some(alternative_names) = self.alternative_names() {
            write!(f, "\t{}:{}", Tag::AlternativeNames, alternative_names)?;
        }

        if let Some(assembly_id) = self.assembly_id() {
            write!(f, "\t{}:{}", Tag::AssemblyId, assembly_id)?;
        }

        if let Some(description) = self.description() {
            write!(f, "\t{}:{}", Tag::Description, description)?;
        }

        if let Some(md5_checksum) = self.md5_checksum() {
            write!(f, "\t{}:{}", Tag::Md5Checksum, md5_checksum)?;
        }

        if let Some(species) = self.species() {
            write!(f, "\t{}:{}", Tag::Species, species)?;
        }

        if let Some(molecule_topology) = self.molecule_topology() {
            write!(f, "\t{}:{}", Tag::MoleculeTopology, molecule_topology)?;
        }

        if let Some(uri) = self.uri() {
            write!(f, "\t{}:{}", Tag::Uri, uri)?;
        }

        for (tag, value) in &self.fields {
            write!(f, "\t{}:{}", tag, value)?;
        }

        Ok(())
    }
}

/// An error returned when a raw SAM header reference sequence fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required tag is missing.
    MissingRequiredTag(Tag),
    /// A tag is invalid.
    InvalidTag(tag::ParseError),
    /// The name tag (`SN`) has an invalid value.
    InvalidName(name::ParseError),
    /// The length tag (`LN`) has an invalid value.
    InvalidLength,
    /// The alternative names tag (`AH`) has an invalid value.
    InvalidAlternativeLocus(alternative_locus::ParseError),
    /// The alternative names tag (`AN`) has an invalid value.
    InvalidAlternativeNames(alternative_names::ParseError),
    /// The MD5 checksum is invalid.
    InvalidMd5Checksum(md5_checksum::ParseError),
    /// The molecule topology is invalid.
    InvalidMoleculeTopology(molecule_topology::ParseError),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "invalid tag: {}", e),
            Self::InvalidName(e) => write!(f, "invalid name: {}", e),
            Self::InvalidLength => write!(f, "invalid length"),
            Self::InvalidAlternativeLocus(e) => write!(f, "invalid alternative locus: {}", e),
            Self::InvalidAlternativeNames(e) => write!(f, "invalid alternative names: {}", e),
            Self::InvalidMd5Checksum(e) => write!(f, "invalid MD5 checksum: {}", e),
            Self::InvalidMoleculeTopology(e) => write!(f, "invalid molecule topology: {}", e),
        }
    }
}

impl TryFrom<Record> for ReferenceSequence {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Kind::ReferenceSequence, record::Value::Map(fields)) => parse_map(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_map(raw_fields: Fields) -> Result<ReferenceSequence, TryFromRecordError> {
    use builder::BuildError;

    let mut builder = ReferenceSequence::builder();

    for (raw_tag, value) in raw_fields {
        let tag = raw_tag.parse().map_err(TryFromRecordError::InvalidTag)?;

        builder = match tag {
            Tag::Name => {
                let name = value.parse().map_err(TryFromRecordError::InvalidName)?;
                builder.set_name(name)
            }
            Tag::Length => {
                let len = value
                    .parse()
                    .map_err(|_| TryFromRecordError::InvalidLength)?;

                builder.set_length(len)
            }
            Tag::AlternativeLocus => {
                let alternative_locus = value
                    .parse()
                    .map_err(TryFromRecordError::InvalidAlternativeLocus)?;
                builder.set_alternative_locus(alternative_locus)
            }
            Tag::AlternativeNames => {
                let alternative_names = value
                    .parse()
                    .map_err(TryFromRecordError::InvalidAlternativeNames)?;
                builder.set_alternative_names(alternative_names)
            }
            Tag::AssemblyId => builder.set_assembly_id(value),
            Tag::Description => builder.set_description(value),
            Tag::Md5Checksum => {
                let md5_checksum = value
                    .parse()
                    .map_err(TryFromRecordError::InvalidMd5Checksum)?;
                builder.set_md5_checksum(md5_checksum)
            }
            Tag::Species => builder.set_species(value),
            Tag::MoleculeTopology => {
                let molecule_topology = value
                    .parse()
                    .map_err(TryFromRecordError::InvalidMoleculeTopology)?;
                builder.set_molecule_topology(molecule_topology)
            }
            Tag::Uri => builder.set_uri(value),
            Tag::Other(..) => builder.insert(tag, value),
        }
    }

    match builder.build() {
        Ok(rs) => Ok(rs),
        Err(BuildError::MissingName) => Err(TryFromRecordError::MissingRequiredTag(Tag::Name)),
        Err(BuildError::MissingLength) => Err(TryFromRecordError::MissingRequiredTag(Tag::Length)),
        Err(BuildError::InvalidLength(_)) => Err(TryFromRecordError::InvalidLength),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() -> Result<(), name::ParseError> {
        assert_eq!(
            ReferenceSequence::new("sq0".parse()?, 0),
            Err(NewError::InvalidLength(0))
        );

        Ok(())
    }

    #[test]
    fn test_fmt() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = ReferenceSequence::builder()
            .set_name("sq0".parse()?)
            .set_length(13)
            .set_md5_checksum(Md5Checksum::from([
                0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6,
                0x15, 0x34,
            ]))
            .build()?;

        assert_eq!(
            reference_sequence.to_string(),
            "@SQ\tSN:sq0\tLN:13\tM5:d7eba311421bbc9d3ada44709dd61534"
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_reference_sequence_with_invalid_record() {
        let record = Record::new(
            record::Kind::Comment,
            record::Value::String(String::from("noodles-sam")),
        );

        assert_eq!(
            ReferenceSequence::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_reference_sequence_with_missing_name(
    ) -> Result<(), record::value::TryFromIteratorError> {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::try_from_iter([
                ("LN", "1"),
                ("M5", "d7eba311421bbc9d3ada44709dd61534"),
            ])?,
        );

        assert_eq!(
            ReferenceSequence::try_from(record),
            Err(TryFromRecordError::MissingRequiredTag(Tag::Name))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_reference_sequence_with_missing_length(
    ) -> Result<(), record::value::TryFromIteratorError> {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::try_from_iter([
                ("SN", "sq0"),
                ("M5", "d7eba311421bbc9d3ada44709dd61534"),
            ])?,
        );

        assert_eq!(
            ReferenceSequence::try_from(record),
            Err(TryFromRecordError::MissingRequiredTag(Tag::Length))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_reference_sequence_with_missing_name_and_length(
    ) -> Result<(), record::value::TryFromIteratorError> {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::try_from_iter([("M5", "d7eba311421bbc9d3ada44709dd61534")])?,
        );

        assert_eq!(
            ReferenceSequence::try_from(record),
            Err(TryFromRecordError::MissingRequiredTag(Tag::Name))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_reference_sequence_with_invalid_name(
    ) -> Result<(), record::value::TryFromIteratorError> {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::try_from_iter([("SN", "*")])?,
        );

        assert!(matches!(
            ReferenceSequence::try_from(record),
            Err(TryFromRecordError::InvalidName(_))
        ));

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_reference_sequence_with_invalid_length(
    ) -> Result<(), record::value::TryFromIteratorError> {
        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::try_from_iter([("LN", "thirteen")])?,
        );

        assert!(matches!(
            ReferenceSequence::try_from(record),
            Err(TryFromRecordError::InvalidLength)
        ));

        let record = Record::new(
            record::Kind::ReferenceSequence,
            record::Value::try_from_iter([("SN", "sq0"), ("LN", "0")])?,
        );

        assert!(matches!(
            ReferenceSequence::try_from(record),
            Err(TryFromRecordError::InvalidLength)
        ));

        Ok(())
    }
}
