//! CRAM record feature.

pub mod code;
pub mod substitution;

pub use self::code::Code;

use noodles_core::Position;

/// A CRAM record feature.
#[allow(missing_docs)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Feature {
    /// A stretch of bases.
    Bases { position: Position, bases: Vec<u8> },
    /// A stretch of quality scores.
    Scores {
        position: Position,
        quality_scores: Vec<u8>,
    },
    /// A base-quality score pair.
    ReadBase {
        position: Position,
        base: u8,
        quality_score: u8,
    },
    /// A base substitution.
    Substitution {
        position: Position,
        value: substitution::Value,
    },
    /// Inserted bases.
    Insertion { position: Position, bases: Vec<u8> },
    /// A number of deleted bases.
    Deletion { position: Position, len: usize },
    /// A single inserted base.
    InsertBase { position: Position, base: u8 },
    /// A single quality score.
    QualityScore {
        position: Position,
        quality_score: u8,
    },
    /// A number of skipped bases.
    ReferenceSkip { position: Position, len: usize },
    /// Soft clipped bases.
    SoftClip { position: Position, bases: Vec<u8> },
    /// A number of padded bases.
    Padding { position: Position, len: usize },
    /// A number of hard clipped bases.
    HardClip { position: Position, len: usize },
}

impl Feature {
    /// Returns the feature code.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::record::{feature::Code, Feature};
    ///
    /// let position = Position::try_from(8)?;
    /// let feature = Feature::Padding{ position, len: 13 };
    ///
    /// assert_eq!(feature.code(), Code::Padding);
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn code(&self) -> Code {
        match self {
            Self::Bases { .. } => Code::Bases,
            Self::Scores { .. } => Code::Scores,
            Self::ReadBase { .. } => Code::ReadBase,
            Self::Substitution { .. } => Code::Substitution,
            Self::Insertion { .. } => Code::Insertion,
            Self::Deletion { .. } => Code::Deletion,
            Self::InsertBase { .. } => Code::InsertBase,
            Self::QualityScore { .. } => Code::QualityScore,
            Self::ReferenceSkip { .. } => Code::ReferenceSkip,
            Self::SoftClip { .. } => Code::SoftClip,
            Self::Padding { .. } => Code::Padding,
            Self::HardClip { .. } => Code::HardClip,
        }
    }

    /// Returns the feature position.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_cram::record::Feature;
    ///
    /// let position = Position::try_from(8)?;
    /// let feature = Feature::Padding{ position, len: 13 };
    ///
    /// assert_eq!(feature.position(), position);
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn position(&self) -> Position {
        match self {
            Self::Bases { position, .. } => *position,
            Self::Scores { position, .. } => *position,
            Self::ReadBase { position, .. } => *position,
            Self::Substitution { position, .. } => *position,
            Self::Insertion { position, .. } => *position,
            Self::Deletion { position, .. } => *position,
            Self::InsertBase { position, .. } => *position,
            Self::QualityScore { position, .. } => *position,
            Self::ReferenceSkip { position, .. } => *position,
            Self::SoftClip { position, .. } => *position,
            Self::Padding { position, .. } => *position,
            Self::HardClip { position, .. } => *position,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_code() {
        let position = Position::MIN;

        assert_eq!(
            Feature::Bases {
                position,
                bases: Vec::new()
            }
            .code(),
            Code::Bases
        );
        assert_eq!(
            Feature::Scores {
                position,
                quality_scores: Vec::new()
            }
            .code(),
            Code::Scores
        );
        assert_eq!(
            Feature::ReadBase {
                position,
                base: b'N',
                quality_score: 0
            }
            .code(),
            Code::ReadBase
        );
        assert_eq!(
            Feature::Substitution {
                position,
                value: substitution::Value::Code(0)
            }
            .code(),
            Code::Substitution
        );
        assert_eq!(
            Feature::Insertion {
                position,
                bases: Vec::new()
            }
            .code(),
            Code::Insertion
        );
        assert_eq!(
            Feature::Deletion { position, len: 0 }.code(),
            Code::Deletion
        );
        assert_eq!(
            Feature::InsertBase {
                position,
                base: b'N'
            }
            .code(),
            Code::InsertBase
        );
        assert_eq!(
            Feature::QualityScore {
                position,
                quality_score: 0
            }
            .code(),
            Code::QualityScore
        );
        assert_eq!(
            Feature::ReferenceSkip { position, len: 0 }.code(),
            Code::ReferenceSkip
        );
        assert_eq!(
            Feature::SoftClip {
                position,
                bases: Vec::new()
            }
            .code(),
            Code::SoftClip
        );
        assert_eq!(Feature::Padding { position, len: 0 }.code(), Code::Padding);
        assert_eq!(
            Feature::HardClip { position, len: 0 }.code(),
            Code::HardClip
        );
    }

    #[test]
    fn test_position() {
        let position = Position::MIN;

        assert_eq!(
            Feature::Bases {
                position,
                bases: Vec::new()
            }
            .position(),
            position
        );
        assert_eq!(
            Feature::Scores {
                position,
                quality_scores: Vec::new()
            }
            .position(),
            position
        );
        assert_eq!(
            Feature::ReadBase {
                position,
                base: b'N',
                quality_score: 0
            }
            .position(),
            position
        );
        assert_eq!(
            Feature::Substitution {
                position,
                value: substitution::Value::Code(0)
            }
            .position(),
            position
        );
        assert_eq!(
            Feature::Insertion {
                position,
                bases: Vec::new()
            }
            .position(),
            position
        );
        assert_eq!(Feature::Deletion { position, len: 0 }.position(), position);
        assert_eq!(
            Feature::InsertBase {
                position,
                base: b'N'
            }
            .position(),
            position
        );
        assert_eq!(
            Feature::QualityScore {
                position,
                quality_score: 0
            }
            .position(),
            position
        );
        assert_eq!(
            Feature::ReferenceSkip { position, len: 0 }.position(),
            position
        );
        assert_eq!(
            Feature::SoftClip {
                position,
                bases: Vec::new()
            }
            .position(),
            position
        );
        assert_eq!(Feature::Padding { position, len: 0 }.position(), position);
        assert_eq!(Feature::HardClip { position, len: 0 }.position(), position);
    }
}
