//! BED record and fields.

pub mod builder;
pub mod color;
pub mod name;
pub mod score;
pub mod strand;

pub use self::{builder::Builder, color::Color, name::Name, score::Score, strand::Strand};

use std::{
    fmt::{self, Write},
    ops::Deref,
};

use noodles_core::Position;

const DELIMITER: char = '\t';

type Block = (usize, usize);

#[derive(Clone, Debug, Eq, PartialEq)]
struct StandardFields {
    reference_sequence_name: String,
    start_position: Position,
    end_position: Position,
    name: Option<Name>,
    score: Option<Score>,
    strand: Option<Strand>,
    thick_start: Position,
    thick_end: Position,
    color: Option<Color>,
    blocks: Vec<Block>,
}

impl StandardFields {
    fn new<N>(reference_sequence_name: N, start_position: Position, end_position: Position) -> Self
    where
        N: Into<String>,
    {
        Self {
            reference_sequence_name: reference_sequence_name.into(),
            start_position,
            end_position,
            name: None,
            score: None,
            strand: None,
            thick_start: start_position,
            thick_end: end_position,
            color: None,
            blocks: Vec::new(),
        }
    }
}

/// Raw BED record optional fields.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct OptionalFields(Vec<String>);

impl Deref for OptionalFields {
    type Target = [String];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for OptionalFields {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, field) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_char(DELIMITER)?;
            }

            f.write_str(field)?;
        }

        Ok(())
    }
}

impl From<Vec<String>> for OptionalFields {
    fn from(fields: Vec<String>) -> Self {
        Self(fields)
    }
}

#[cfg(test)]
mod optional_fields_tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let fields = OptionalFields::default();
        assert_eq!(fields.to_string(), "");

        let fields = OptionalFields::from(vec![String::from("n")]);
        assert_eq!(fields.to_string(), "n");

        let fields = OptionalFields::from(vec![String::from("n"), String::from("d")]);
        assert_eq!(fields.to_string(), "n\td");
    }
}

/// A BED record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RecordBuf<const N: u8> {
    standard_fields: StandardFields,
    optional_fields: OptionalFields,
}

/// A trait describing the number of standard fields in a BED record.
pub trait BedN<const N: u8> {}

impl BedN<3> for RecordBuf<3> {}
impl BedN<3> for RecordBuf<4> {}
impl BedN<3> for RecordBuf<5> {}
impl BedN<3> for RecordBuf<6> {}
impl BedN<3> for RecordBuf<7> {}
impl BedN<3> for RecordBuf<8> {}
impl BedN<3> for RecordBuf<9> {}
impl BedN<3> for RecordBuf<12> {}

impl BedN<4> for RecordBuf<4> {}
impl BedN<4> for RecordBuf<5> {}
impl BedN<4> for RecordBuf<6> {}
impl BedN<4> for RecordBuf<7> {}
impl BedN<4> for RecordBuf<8> {}
impl BedN<4> for RecordBuf<9> {}
impl BedN<4> for RecordBuf<12> {}

impl BedN<5> for RecordBuf<5> {}
impl BedN<5> for RecordBuf<6> {}
impl BedN<5> for RecordBuf<7> {}
impl BedN<5> for RecordBuf<8> {}
impl BedN<5> for RecordBuf<9> {}
impl BedN<5> for RecordBuf<12> {}

impl BedN<6> for RecordBuf<6> {}
impl BedN<6> for RecordBuf<7> {}
impl BedN<6> for RecordBuf<8> {}
impl BedN<6> for RecordBuf<9> {}
impl BedN<6> for RecordBuf<12> {}

impl BedN<7> for RecordBuf<7> {}
impl BedN<7> for RecordBuf<8> {}
impl BedN<7> for RecordBuf<9> {}
impl BedN<7> for RecordBuf<12> {}

impl BedN<8> for RecordBuf<8> {}
impl BedN<8> for RecordBuf<9> {}
impl BedN<8> for RecordBuf<12> {}

impl BedN<9> for RecordBuf<9> {}
impl BedN<9> for RecordBuf<12> {}

impl BedN<12> for RecordBuf<12> {}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<3>,
{
    /// Creates a BED record builder.
    pub fn builder() -> Builder<N> {
        Builder::default()
    }

    fn new(standard_fields: StandardFields, optional_fields: OptionalFields) -> Self {
        Self {
            standard_fields,
            optional_fields,
        }
    }

    /// Returns the reference sequence name (`chrom`).
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
    pub fn reference_sequence_name(&self) -> &str {
        &self.standard_fields.reference_sequence_name
    }

    /// Returns the feature start position (`chromStart`).
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
    pub fn start_position(&self) -> Position {
        self.standard_fields.start_position
    }

    /// Returns the feature end position (`chromEnd`).
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
    pub fn end_position(&self) -> Position {
        self.standard_fields.end_position
    }

    /// Returns the list of raw optional fields.
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
    /// assert!(record.optional_fields().is_empty());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn optional_fields(&self) -> &OptionalFields {
        &self.optional_fields
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<4>,
{
    /// Returns the feature name (`name`).
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
    pub fn name(&self) -> Option<&Name> {
        self.standard_fields.name.as_ref()
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<5>,
{
    /// Returns the score (`score`).
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
    pub fn score(&self) -> Option<Score> {
        self.standard_fields.score
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<6>,
{
    /// Returns the feature strand (`strand`).
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
    pub fn strand(&self) -> Option<Strand> {
        self.standard_fields.strand
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<7>,
{
    /// Returns the thick start position (`thickStart`).
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
    pub fn thick_start(&self) -> Position {
        self.standard_fields.thick_start
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<8>,
{
    /// Returns the thick end position (`thickEnd`).
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
    pub fn thick_end(&self) -> Position {
        self.standard_fields.thick_end
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<9>,
{
    /// Returns the color (`itemRgb`).
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
    pub fn color(&self) -> Option<Color> {
        self.standard_fields.color
    }
}

impl<const N: u8> RecordBuf<N>
where
    Self: BedN<12>,
{
    /// Returns the blocks (`[(blockStarts, blockSizes)]`).
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
    pub fn blocks(&self) -> &[Block] {
        &self.standard_fields.blocks
    }
}
