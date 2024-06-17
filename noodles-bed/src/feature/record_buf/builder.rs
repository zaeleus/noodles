//! BED record builder.

use std::{error, fmt};

use noodles_core::Position;

use super::{BedN, Block, Color, Name, OptionalFields, RecordBuf, Score, StandardFields, Strand};

/// A BED record builder.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Builder<const N: u8> {
    reference_sequence_name: Option<String>,
    start_position: Option<Position>,
    end_position: Option<Position>,
    name: Option<Name>,
    score: Option<Score>,
    strand: Option<Strand>,
    thick_start: Option<Position>,
    thick_end: Option<Position>,
    color: Option<Color>,
    blocks: Vec<Block>,
    optional_fields: OptionalFields,
}

impl BedN<3> for Builder<3> {}
impl BedN<3> for Builder<4> {}
impl BedN<3> for Builder<5> {}
impl BedN<3> for Builder<6> {}
impl BedN<3> for Builder<7> {}
impl BedN<3> for Builder<8> {}
impl BedN<3> for Builder<9> {}
impl BedN<3> for Builder<12> {}

impl BedN<4> for Builder<4> {}
impl BedN<4> for Builder<5> {}
impl BedN<4> for Builder<6> {}
impl BedN<4> for Builder<7> {}
impl BedN<4> for Builder<8> {}
impl BedN<4> for Builder<9> {}
impl BedN<4> for Builder<12> {}

impl BedN<5> for Builder<5> {}
impl BedN<5> for Builder<6> {}
impl BedN<5> for Builder<7> {}
impl BedN<5> for Builder<8> {}
impl BedN<5> for Builder<9> {}
impl BedN<5> for Builder<12> {}

impl BedN<6> for Builder<6> {}
impl BedN<6> for Builder<7> {}
impl BedN<6> for Builder<8> {}
impl BedN<6> for Builder<9> {}
impl BedN<6> for Builder<12> {}

impl BedN<7> for Builder<7> {}
impl BedN<7> for Builder<8> {}
impl BedN<7> for Builder<9> {}
impl BedN<7> for Builder<12> {}

impl BedN<8> for Builder<8> {}
impl BedN<8> for Builder<9> {}
impl BedN<8> for Builder<12> {}

impl BedN<9> for Builder<9> {}
impl BedN<9> for Builder<12> {}

impl BedN<12> for Builder<12> {}

impl<const N: u8> Builder<N>
where
    Self: BedN<3>,
{
    /// Sets the reference sequence name (`chrom`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.reference_sequence_name(), "sq0");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_reference_sequence_name<M>(mut self, reference_sequence_name: M) -> Self
    where
        M: Into<String>,
    {
        self.reference_sequence_name = Some(reference_sequence_name.into());
        self
    }

    /// Sets the feature start position (`chromStart`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let start_position = Position::try_from(8)?;
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(start_position)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.start_position(), start_position);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_start_position(mut self, start_position: Position) -> Self {
        self.start_position = Some(start_position);
        self
    }

    /// Sets the feature end position (`chromEnd`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let end_position = Position::try_from(13)?;
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(end_position)
    ///     .build()?;
    ///
    /// assert_eq!(record.end_position(), end_position);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_end_position(mut self, end_position: Position) -> Self {
        self.end_position = Some(end_position);
        self
    }

    /// Sets the list of raw optional fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::OptionalFields};
    /// use noodles_core::Position;
    ///
    /// let optional_fields = OptionalFields::from(vec![String::from("n")]);
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_optional_fields(optional_fields.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.optional_fields(), &optional_fields);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_optional_fields(mut self, optional_fields: OptionalFields) -> Self {
        self.optional_fields = optional_fields;
        self
    }
}

impl Builder<3> {
    /// Builds a BED3 record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<3>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<RecordBuf<3>, BuildError> {
        let reference_sequence_name = self
            .reference_sequence_name
            .ok_or(BuildError::MissingReferenceSequenceName)?;

        let start_position = self
            .start_position
            .ok_or(BuildError::MissingStartPosition)?;

        let end_position = self.end_position.ok_or(BuildError::MissingEndPosition)?;

        let standard_fields =
            StandardFields::new(reference_sequence_name, start_position, end_position);

        Ok(RecordBuf::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<4>,
{
    /// Sets the feature name (`name`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::Name};
    /// use noodles_core::Position;
    ///
    /// let name = Name::from("ndls1");
    ///
    /// let record = bed::feature::RecordBuf::<4>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_name(name.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.name(), Some(&name));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_name(mut self, name: Name) -> Self {
        self.name = Some(name);
        self
    }
}

impl Builder<4> {
    /// Builds a BED4 record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<4>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<RecordBuf<4>, BuildError> {
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

        Ok(RecordBuf::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<5>,
{
    /// Sets the score (`score`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::Score};
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<5>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_score(Score::try_from(21)?)
    ///     .build()?;
    ///
    /// assert_eq!(record.score().map(u16::from), Some(21));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_score(mut self, score: Score) -> Self {
        self.score = Some(score);
        self
    }
}

impl Builder<5> {
    /// Builds a BED5 record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<5>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<RecordBuf<5>, BuildError> {
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

        Ok(RecordBuf::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<6>,
{
    /// Sets the feature strand (`strand`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::Strand};
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<6>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_strand(Strand::Forward)
    ///     .build()?;
    ///
    /// assert_eq!(record.strand(), Some(Strand::Forward));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }
}

impl Builder<6> {
    /// Builds a BED6 record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<6>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<RecordBuf<6>, BuildError> {
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

        Ok(RecordBuf::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<7>,
{
    /// Sets the thick start position (`thickStart`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let thick_start = Position::try_from(8)?;
    ///
    /// let record = bed::feature::RecordBuf::<7>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_thick_start(thick_start)
    ///     .build()?;
    ///
    /// assert_eq!(record.thick_start(), thick_start);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_thick_start(mut self, thick_start: Position) -> Self {
        self.thick_start = Some(thick_start);
        self
    }
}

impl Builder<7> {
    /// Builds a BED7 record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<7>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<RecordBuf<7>, BuildError> {
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
        standard_fields.thick_start = self.thick_start.unwrap_or(start_position);

        Ok(RecordBuf::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<8>,
{
    /// Sets the thick end position (`thickEnd`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let thick_end = Position::try_from(13)?;
    ///
    /// let record = bed::feature::RecordBuf::<8>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_thick_end(thick_end)
    ///     .build()?;
    ///
    /// assert_eq!(record.thick_end(), thick_end);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_thick_end(mut self, thick_end: Position) -> Self {
        self.thick_end = Some(thick_end);
        self
    }
}

impl Builder<8> {
    /// Builds a BED8 record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<8>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<RecordBuf<8>, BuildError> {
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
        standard_fields.thick_start = self.thick_start.unwrap_or(start_position);
        standard_fields.thick_end = self.thick_end.unwrap_or(end_position);

        Ok(RecordBuf::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<9>,
{
    /// Sets the color (`itemRgb`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed::{self as bed, feature::record_buf::Color};
    /// use noodles_core::Position;
    ///
    /// let thick_end = Position::try_from(13)?;
    ///
    /// let record = bed::feature::RecordBuf::<9>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_color(Color::RED)
    ///     .build()?;
    ///
    /// assert_eq!(record.color(), Some(Color::RED));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_color(mut self, color: Color) -> Self {
        self.color = Some(color);
        self
    }
}

impl Builder<9> {
    /// Builds a BED9 record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<9>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<RecordBuf<9>, BuildError> {
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
        standard_fields.thick_start = self.thick_start.unwrap_or(start_position);
        standard_fields.thick_end = self.thick_end.unwrap_or(end_position);
        standard_fields.color = self.color;

        Ok(RecordBuf::new(standard_fields, self.optional_fields))
    }
}

impl<const N: u8> Builder<N>
where
    Self: BedN<12>,
{
    /// Sets the blocks (`[(blockStarts, blockSizes)]`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let blocks = vec![(0, 2)];
    ///
    /// let record = bed::feature::RecordBuf::<12>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .set_blocks(blocks.clone())
    ///     .build()?;
    ///
    /// assert_eq!(record.blocks(), &blocks);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_blocks(mut self, blocks: Vec<Block>) -> Self {
        self.blocks = blocks;
        self
    }
}

impl Builder<12> {
    /// Builds a BED12 record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// use noodles_core::Position;
    ///
    /// let record = bed::feature::RecordBuf::<12>::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_start_position(Position::try_from(8)?)
    ///     .set_end_position(Position::try_from(13)?)
    ///     .build()?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn build(self) -> Result<RecordBuf<12>, BuildError> {
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
        standard_fields.thick_start = self.thick_start.unwrap_or(start_position);
        standard_fields.thick_end = self.thick_end.unwrap_or(end_position);
        standard_fields.color = self.color;
        standard_fields.blocks = self.blocks;

        Ok(RecordBuf::new(standard_fields, self.optional_fields))
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
