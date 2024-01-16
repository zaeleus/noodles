//! SAM header reference sequence builder.

use std::num::NonZeroUsize;

use super::{AlternativeLocus, AlternativeNames, Md5Checksum, MoleculeTopology, ReferenceSequence};
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header reference sequence builder.
#[derive(Debug, Default)]
pub struct Builder {
    length: Option<NonZeroUsize>,
    alternative_locus: Option<AlternativeLocus>,
    alternative_names: Option<AlternativeNames>,
    assembly_id: Option<Vec<u8>>,
    description: Option<Vec<u8>>,
    md5_checksum: Option<Md5Checksum>,
    species: Option<Vec<u8>>,
    molecule_topology: Option<MoleculeTopology>,
    uri: Option<Vec<u8>>,
}

impl map::Builder<ReferenceSequence> {
    /// Sets a reference sequence length.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let length = NonZeroUsize::try_from(13)?;
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(length)
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.length(), length);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_length(mut self, length: NonZeroUsize) -> Self {
        self.inner.length = Some(length);
        self
    }

    /// Sets an alternative locus.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::AlternativeLocus, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(NonZeroUsize::try_from(13)?)
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
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::AlternativeNames, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let alternative_names: AlternativeNames = "0,SQ.0".parse()?;
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(NonZeroUsize::try_from(13)?)
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
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(NonZeroUsize::try_from(13)?)
    ///     .set_assembly_id("ref")
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.assembly_id(), Some(&b"ref"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_assembly_id<I>(mut self, assembly_id: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.assembly_id = Some(assembly_id.into());
        self
    }

    /// Sets a description.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(NonZeroUsize::try_from(13)?)
    ///     .set_description("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.description(), Some(&b"noodles"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_description<I>(mut self, description: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.description = Some(description.into());
        self
    }

    /// Sets an MD5 checksum.
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
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(NonZeroUsize::try_from(13)?)
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
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(NonZeroUsize::try_from(13)?)
    ///     .set_species("human")
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.species(), Some(&b"human"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_species<I>(mut self, species: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.species = Some(species.into());
        self
    }

    /// Sets a molecule topology.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{
    ///     map::{reference_sequence::MoleculeTopology, ReferenceSequence},
    ///     Map,
    /// };
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(NonZeroUsize::try_from(13)?)
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
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
    ///
    /// let reference_sequence = Map::<ReferenceSequence>::builder()
    ///     .set_length(NonZeroUsize::try_from(13)?)
    ///     .set_uri("file:///tmp/ref.fasta")
    ///     .build()?;
    ///
    /// assert_eq!(reference_sequence.uri(), Some(&b"file:///tmp/ref.fasta"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_uri<I>(mut self, uri: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.uri = Some(uri.into());
        self
    }
}

impl map::builder::Inner<ReferenceSequence> for Builder {
    fn build(self) -> Result<ReferenceSequence, BuildError> {
        let length = self.length.ok_or(BuildError::MissingField("LN"))?;

        Ok(ReferenceSequence {
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
