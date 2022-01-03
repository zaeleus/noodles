//! BED record builder.

use std::{error, fmt};

use super::{BedN, OptionalFields, Record, Score, StandardFields, Strand};

/// A BED record builder.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Builder<const N: u8> {
    reference_sequence_name: Option<String>,
    start_position: Option<u64>,
    end_position: Option<u64>,
    name: Option<String>,
    score: Option<Score>,
    strand: Option<Strand>,
    optional_fields: OptionalFields,
}

impl BedN<3> for Builder<3> {}
impl BedN<3> for Builder<4> {}
impl BedN<3> for Builder<5> {}
impl BedN<3> for Builder<6> {}

impl BedN<4> for Builder<4> {}
impl BedN<4> for Builder<5> {}
impl BedN<4> for Builder<6> {}

impl BedN<5> for Builder<5> {}
impl BedN<5> for Builder<6> {}

impl BedN<6> for Builder<6> {}

impl<const N: u8> Builder<N>
where
    Self: BedN<3>,
{
    /// Sets the reference sequence name (`chrom`).
    pub fn set_reference_sequence_name<M>(mut self, reference_sequence_name: M) -> Self
    where
        M: Into<String>,
    {
        self.reference_sequence_name = Some(reference_sequence_name.into());
        self
    }

    /// Sets the feature start position (`chromStart`).
    pub fn set_start_position(mut self, start_position: u64) -> Self {
        self.start_position = Some(start_position);
        self
    }

    /// Sets the feature end position (`chromEnd`).
    pub fn set_end_position(mut self, end_position: u64) -> Self {
        self.end_position = Some(end_position);
        self
    }

    /// Sets the list of raw optional fields.
    pub fn set_optional_fields(mut self, optional_fields: OptionalFields) -> Self {
        self.optional_fields = optional_fields;
        self
    }
}

impl Builder<3> {
    /// Builds a BED3 record.
    pub fn build(self) -> Result<Record<3>, BuildError> {
        let reference_sequence_name = self
            .reference_sequence_name
            .ok_or(BuildError::MissingReferenceSequenceName)?;

        let start_position = self
            .start_position
            .ok_or(BuildError::MissingStartPosition)?;

        let end_position = self.end_position.ok_or(BuildError::MissingEndPosition)?;

        let standard_fields =
            StandardFields::new(reference_sequence_name, start_position, end_position);

        Ok(Record::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<4>,
{
    /// Sets the feature name (`name`).
    pub fn set_name<M>(mut self, name: M) -> Self
    where
        M: Into<String>,
    {
        self.name = Some(name.into());
        self
    }
}

impl Builder<4> {
    /// Builds a BED4 record.
    pub fn build(self) -> Result<Record<4>, BuildError> {
        let reference_sequence_name = self
            .reference_sequence_name
            .ok_or(BuildError::MissingReferenceSequenceName)?;

        let start_position = self
            .start_position
            .ok_or(BuildError::MissingStartPosition)?;

        let end_position = self.end_position.ok_or(BuildError::MissingEndPosition)?;

        let mut standard_fields =
            StandardFields::new(reference_sequence_name, start_position, end_position);
        standard_fields.name = self.name;

        Ok(Record::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<5>,
{
    /// Sets the score (`score`).
    pub fn set_score(mut self, score: Score) -> Self {
        self.score = Some(score);
        self
    }
}

impl Builder<5> {
    /// Builds a BED5 record.
    pub fn build(self) -> Result<Record<5>, BuildError> {
        let reference_sequence_name = self
            .reference_sequence_name
            .ok_or(BuildError::MissingReferenceSequenceName)?;

        let start_position = self
            .start_position
            .ok_or(BuildError::MissingStartPosition)?;

        let end_position = self.end_position.ok_or(BuildError::MissingEndPosition)?;

        let mut standard_fields =
            StandardFields::new(reference_sequence_name, start_position, end_position);
        standard_fields.name = self.name;
        standard_fields.score = self.score;

        Ok(Record::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<6>,
{
    /// Sets the feature strand (`strand`).
    pub fn set_strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }
}

impl Builder<6> {
    /// Builds a BED6 record.
    pub fn build(self) -> Result<Record<6>, BuildError> {
        let reference_sequence_name = self
            .reference_sequence_name
            .ok_or(BuildError::MissingReferenceSequenceName)?;

        let start_position = self
            .start_position
            .ok_or(BuildError::MissingStartPosition)?;

        let end_position = self.end_position.ok_or(BuildError::MissingEndPosition)?;

        let mut standard_fields =
            StandardFields::new(reference_sequence_name, start_position, end_position);
        standard_fields.name = self.name;
        standard_fields.score = self.score;
        standard_fields.strand = self.strand;

        Ok(Record::new(standard_fields, self.optional_fields))
    }
}

/// An error returned when a BED record fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
    /// The start position is missing.
    MissingStartPosition,
    /// The end position is missing.
    MissingEndPosition,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceName => f.write_str("missing reference sequence name"),
            Self::MissingStartPosition => f.write_str("missing start position"),
            Self::MissingEndPosition => f.write_str("missing end position"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_3() -> Result<(), BuildError> {
        Builder::<3>::default()
            .set_reference_sequence_name("sq0")
            .set_start_position(8)
            .set_end_position(13)
            .build()?;

        Ok(())
    }

    #[test]
    fn test_build_4() -> Result<(), BuildError> {
        Builder::<4>::default()
            .set_reference_sequence_name("sq0")
            .set_start_position(8)
            .set_end_position(13)
            .set_name("ndls1")
            .build()?;

        Ok(())
    }

    #[test]
    fn test_build_5() -> Result<(), Box<dyn std::error::Error>> {
        Builder::<5>::default()
            .set_reference_sequence_name("sq0")
            .set_start_position(8)
            .set_end_position(13)
            .set_score(Score::try_from(21)?)
            .build()?;

        Ok(())
    }

    #[test]
    fn test_build_6() -> Result<(), BuildError> {
        Builder::<6>::default()
            .set_reference_sequence_name("sq0")
            .set_start_position(8)
            .set_end_position(13)
            .set_strand(Strand::Forward)
            .build()?;

        Ok(())
    }
}
