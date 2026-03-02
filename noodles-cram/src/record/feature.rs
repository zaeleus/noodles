pub mod code;

pub use self::code::Code;

use std::borrow::Cow;

use noodles_core::Position;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Feature<'c> {
    Bases {
        position: Position,
        bases: Cow<'c, [u8]>,
    },
    Scores {
        position: Position,
        quality_scores: Cow<'c, [u8]>,
    },
    ReadBase {
        position: Position,
        base: u8,
        quality_score: u8,
    },
    Substitution {
        position: Position,
        code: u8,
    },
    Insertion {
        position: Position,
        bases: Cow<'c, [u8]>,
    },
    Deletion {
        position: Position,
        len: usize,
    },
    InsertBase {
        position: Position,
        base: u8,
    },
    QualityScore {
        position: Position,
        quality_score: u8,
    },
    ReferenceSkip {
        position: Position,
        len: usize,
    },
    SoftClip {
        position: Position,
        bases: Cow<'c, [u8]>,
    },
    Padding {
        position: Position,
        len: usize,
    },
    HardClip {
        position: Position,
        len: usize,
    },
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
                bases: Cow::Borrowed(&[])
            }
            .position(),
            position
        );
        assert_eq!(
            Feature::Scores {
                position,
                quality_scores: Cow::Borrowed(&[])
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
                bases: Cow::Borrowed(&[])
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
                bases: Cow::Borrowed(&[])
            }
            .position(),
            position
        );
        assert_eq!(Feature::Padding { position, len: 0 }.position(), position);
        assert_eq!(Feature::HardClip { position, len: 0 }.position(), position);
    }
}
