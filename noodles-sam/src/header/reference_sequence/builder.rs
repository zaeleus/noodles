//! SAM header reference sequence builder.

use std::{collections::HashMap, error, fmt, num::NonZeroUsize};

use super::{
    AlternativeLocus, AlternativeNames, Md5Checksum, MoleculeTopology, Name, ReferenceSequence, Tag,
};

/// A SAM header reference sequence builder.
#[derive(Debug, Default)]
pub struct Builder {
    name: Option<Name>,
    len: Option<usize>,
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

/// An error returned when a SAM header reference sequence fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The name is missing.
    MissingName,
    /// The length is missing.
    MissingLength,
    /// The length is invalid.
    InvalidLength(usize),
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingName => f.write_str("missing name"),
            Self::MissingLength => f.write_str("missing length"),
            Self::InvalidLength(len) => write!(f, "invalid length: {}", len),
        }
    }
}

impl Builder {
    /// Sets a reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .build()?;
    ///
    /// assert_eq!(**reference_sequence.name(), "sq0");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_name(mut self, name: Name) -> Self {
        self.name = Some(name);
        self
    }

    /// Sets a reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .build()?;
    ///
    /// assert_eq!(usize::from(reference_sequence.len()), 13);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_length(mut self, len: usize) -> Self {
        self.len = Some(len);
        self
    }

    /// Sets an alternative locus.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::AlternativeLocus, ReferenceSequence};
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_alternative_locus(AlternativeLocus::Unknown)
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     reference_sequence.alternative_locus(),
    ///     Some(&AlternativeLocus::Unknown)
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_alternative_locus(mut self, alternative_locus: AlternativeLocus) -> Self {
        self.alternative_locus = Some(alternative_locus);
        self
    }

    /// Sets alternative names.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{
    ///     reference_sequence::AlternativeNames,
    ///     ReferenceSequence,
    /// };
    ///
    /// let alternative_names: AlternativeNames = "0,SQ.0".parse()?;
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_alternative_names(alternative_names.clone())
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.alternative_names(), Some(&alternative_names));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_alternative_names(mut self, alternative_names: AlternativeNames) -> Self {
        self.alternative_names = Some(alternative_names);
        self
    }

    /// Sets a genome assembly ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_assembly_id("ref")
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.assembly_id(), Some("ref"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_assembly_id<I>(mut self, assembly_id: I) -> Self
    where
        I: Into<String>,
    {
        self.assembly_id = Some(assembly_id.into());
        self
    }

    /// Sets a description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_description("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.description(), Some("noodles"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_description<I>(mut self, description: I) -> Self
    where
        I: Into<String>,
    {
        self.description = Some(description.into());
        self
    }

    /// Sets an MD5 checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::Md5Checksum, ReferenceSequence};
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_md5_checksum(Md5Checksum::from([
    ///         0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d,
    ///         0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6, 0x15, 0x34,
    ///     ]))
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.md5_checksum(), Some(Md5Checksum::from([
    ///     0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d,
    ///     0x3a, 0xda, 0x44, 0x70, 0x9d, 0xd6, 0x15, 0x34,
    /// ])));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_md5_checksum(mut self, md5_checksum: Md5Checksum) -> Self {
        self.md5_checksum = Some(md5_checksum);
        self
    }

    /// Sets a species.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_species("human")
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.species(), Some("human"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_species<I>(mut self, species: I) -> Self
    where
        I: Into<String>,
    {
        self.species = Some(species.into());
        self
    }

    /// Sets a molecule topology.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::MoleculeTopology, ReferenceSequence};
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_molecule_topology(MoleculeTopology::Linear)
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.molecule_topology(), Some(MoleculeTopology::Linear));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_molecule_topology(mut self, molecule_topology: MoleculeTopology) -> Self {
        self.molecule_topology = Some(molecule_topology);
        self
    }

    /// Sets a URI.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_uri("file:///tmp/ref.fasta")
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.uri(), Some("file:///tmp/ref.fasta"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_uri<I>(mut self, uri: I) -> Self
    where
        I: Into<String>,
    {
        self.uri = Some(uri.into());
        self
    }

    /// Inserts a tag-raw value pair.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{reference_sequence::Tag, ReferenceSequence};
    ///
    /// let zn = Tag::Other([b'z', b'n']);
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .insert(zn, "noodles")
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     reference_sequence.fields().get(&zn),
    ///     Some(&String::from("noodles"))
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn insert<I>(mut self, tag: Tag, value: I) -> Self
    where
        I: Into<String>,
    {
        self.fields.insert(tag, value.into());
        self
    }

    /// Builds a reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReferenceSequence;
    ///
    /// let reference_sequence = ReferenceSequence::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .build()?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<ReferenceSequence, BuildError> {
        let name = self.name.ok_or(BuildError::MissingName)?;

        let len = self
            .len
            .ok_or(BuildError::MissingLength)
            .and_then(|n| NonZeroUsize::new(n).ok_or(BuildError::InvalidLength(n)))?;

        Ok(ReferenceSequence {
            name,
            len,
            alternative_locus: self.alternative_locus,
            alternative_names: self.alternative_names,
            assembly_id: self.assembly_id,
            description: self.description,
            md5_checksum: self.md5_checksum,
            species: self.species,
            molecule_topology: self.molecule_topology,
            uri: self.uri,
            fields: self.fields,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::reference_sequence::name;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.name.is_none());
        assert!(builder.len.is_none());
        assert!(builder.alternative_locus.is_none());
        assert!(builder.alternative_names.is_none());
        assert!(builder.assembly_id.is_none());
        assert!(builder.description.is_none());
        assert!(builder.md5_checksum.is_none());
        assert!(builder.species.is_none());
        assert!(builder.molecule_topology.is_none());
        assert!(builder.uri.is_none());
        assert!(builder.fields.is_empty());
    }

    #[test]
    fn test_build() -> Result<(), name::ParseError> {
        assert_eq!(
            Builder::default().set_length(13).build(),
            Err(BuildError::MissingName)
        );

        let name: Name = "sq0".parse()?;

        assert_eq!(
            Builder::default().set_name(name.clone()).build(),
            Err(BuildError::MissingLength)
        );

        assert_eq!(
            Builder::default().set_name(name).set_length(0).build(),
            Err(BuildError::InvalidLength(0))
        );

        Ok(())
    }
}
