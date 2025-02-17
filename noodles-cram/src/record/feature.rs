//! CRAM record feature.

pub mod code;

pub use self::code::Code;

use noodles_core::Position;

/// A CRAM record feature.
#[allow(missing_docs)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Feature<'c> {
    /// A stretch of bases.
    Bases { position: Position, bases: &'c [u8] },
    /// A stretch of quality scores.
    Scores {
        position: Position,
        quality_scores: &'c [u8],
    },
    /// A base-quality score pair.
    ReadBase {
        position: Position,
        base: u8,
        quality_score: u8,
    },
    /// A base substitution.
    Substitution { position: Position, code: u8 },
    /// Inserted bases.
    Insertion { position: Position, bases: &'c [u8] },
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
    SoftClip { position: Position, bases: &'c [u8] },
    /// A number of padded bases.
    Padding { position: Position, len: usize },
    /// A number of hard clipped bases.
    HardClip { position: Position, len: usize },
}

impl Feature<'_> {
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
    fn test_position() {
        let position = Position::MIN;

        assert_eq!(
            Feature::Bases {
                position,
                bases: &[]
            }
            .position(),
            position
        );
        assert_eq!(
            Feature::Scores {
                position,
                quality_scores: &[]
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
            Feature::Substitution { position, code: 0 }.position(),
            position
        );
        assert_eq!(
            Feature::Insertion {
                position,
                bases: &[]
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
                bases: &[]
            }
            .position(),
            position
        );
        assert_eq!(Feature::Padding { position, len: 0 }.position(), position);
        assert_eq!(Feature::HardClip { position, len: 0 }.position(), position);
    }
}
