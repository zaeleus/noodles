//! SAM header reference sequence builder.

use std::num::NonZeroUsize;

use super::{
    AlternativeLocus, AlternativeNames, Md5Checksum, MoleculeTopology, Name, ReferenceSequence,
};
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header reference sequence builder.
#[derive(Debug, Default)]
pub struct Builder {
    name: Option<Name>,
    length: Option<usize>,
    alternative_locus: Option<AlternativeLocus>,
    alternative_names: Option<AlternativeNames>,
    assembly_id: Option<String>,
    description: Option<String>,
    md5_checksum: Option<Md5Checksum>,
    species: Option<String>,
    molecule_topology: Option<MoleculeTopology>,
    uri: Option<String>,
}

impl map::Builder<ReferenceSequence> {
    /// Sets a reference sequence name.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .build()?;
    ///
    /// assert_eq!(**reference_sequence.name(), "sq0");
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_name(mut self, name: Name) -> Self {
        self.inner.name = Some(name);
        self
    }

    /// Sets a reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .build()?;
    ///
    /// assert_eq!(usize::from(reference_sequence.length()), 13);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_length(mut self, length: usize) -> Self {
        self.inner.length = Some(length);
        self
    }

    /// Sets an alternative locus.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::AlternativeLocus, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
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
        self.inner.alternative_locus = Some(alternative_locus);
        self
    }

    /// Sets alternative names.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::AlternativeNames, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let alternative_names: AlternativeNames = "0,SQ.0".parse()?;
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_alternative_names(alternative_names.clone())
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.alternative_names(), Some(&alternative_names));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_alternative_names(mut self, alternative_names: AlternativeNames) -> Self {
        self.inner.alternative_names = Some(alternative_names);
        self
    }

    /// Sets a genome assembly ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
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
        self.inner.assembly_id = Some(assembly_id.into());
        self
    }

    /// Sets a description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
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
        self.inner.description = Some(description.into());
        self
    }

    /// Sets an MD5 checksum.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::Md5Checksum, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
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
        self.inner.md5_checksum = Some(md5_checksum);
        self
    }

    /// Sets a species.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
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
        self.inner.species = Some(species.into());
        self
    }

    /// Sets a molecule topology.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::MoleculeTopology, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_name("sq0".parse()?)
    ///     .set_length(13)
    ///     .set_molecule_topology(MoleculeTopology::Linear)
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.molecule_topology(), Some(MoleculeTopology::Linear));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_molecule_topology(mut self, molecule_topology: MoleculeTopology) -> Self {
        self.inner.molecule_topology = Some(molecule_topology);
        self
    }

    /// Sets a URI.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
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
        self.inner.uri = Some(uri.into());
        self
    }
}

impl map::builder::Inner<ReferenceSequence> for Builder {
    fn build(self) -> Result<ReferenceSequence, BuildError> {
        let name = self.name.ok_or(BuildError::MissingField("SN"))?;

        let length = self
            .length
            .ok_or(BuildError::MissingField("LN"))
            .and_then(|n| NonZeroUsize::new(n).ok_or(BuildError::InvalidValue("LN")))?;

        Ok(ReferenceSequence {
            name,
            length,
            alternative_locus: self.alternative_locus,
            alternative_names: self.alternative_names,
            assembly_id: self.assembly_id,
            description: self.description,
            md5_checksum: self.md5_checksum,
            species: self.species,
            molecule_topology: self.molecule_topology,
            uri: self.uri,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.name.is_none());
        assert!(builder.length.is_none());
        assert!(builder.alternative_locus.is_none());
        assert!(builder.alternative_names.is_none());
        assert!(builder.assembly_id.is_none());
        assert!(builder.description.is_none());
        assert!(builder.md5_checksum.is_none());
        assert!(builder.species.is_none());
        assert!(builder.molecule_topology.is_none());
        assert!(builder.uri.is_none());
    }
}
