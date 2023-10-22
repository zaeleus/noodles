//! CRAM record feature.

pub mod code;
pub mod substitution;

pub use self::code::Code;

use noodles_core::Position;
use noodles_sam::record::sequence::Base;

/// A CRAM record feature.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Feature {
    /// A stretch of bases (position, bases).
    Bases(Position, Vec<Base>),
    /// A stretch of quality scores (position, quality scores).
    Scores(Position, Vec<u8>),
    /// A base-quality score pair (position, base, quality score).
    ReadBase(Position, Base, u8),
    /// A base substitution (position, code (read) / base (write)).
    Substitution(Position, substitution::Value),
    /// Inserted bases (position, bases).
    Insertion(Position, Vec<Base>),
    /// A number of deleted bases (position, length).
    Deletion(Position, usize),
    /// A single inserted base (position, base).
    InsertBase(Position, Base),
    /// A single quality score (position, score).
    QualityScore(Position, u8),
    /// A number of skipped bases (position, length).
    ReferenceSkip(Position, usize),
    /// Soft clipped bases (position, bases).
    SoftClip(Position, Vec<Base>),
    /// A number of padded bases (position, length).
    Padding(Position, usize),
    /// A number of hard clipped bases (position, length).
    HardClip(Position, usize),
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
    /// let feature = Feature::Padding(position, 13);
    ///
    /// assert_eq!(feature.code(), Code::Padding);
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn code(&self) -> Code {
        match self {
            Self::Bases(..) => Code::Bases,
            Self::Scores(..) => Code::Scores,
            Self::ReadBase(..) => Code::ReadBase,
            Self::Substitution(..) => Code::Substitution,
            Self::Insertion(..) => Code::Insertion,
            Self::Deletion(..) => Code::Deletion,
            Self::InsertBase(..) => Code::InsertBase,
            Self::QualityScore(..) => Code::QualityScore,
            Self::ReferenceSkip(..) => Code::ReferenceSkip,
            Self::SoftClip(..) => Code::SoftClip,
            Self::Padding(..) => Code::Padding,
            Self::HardClip(..) => Code::HardClip,
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
    /// let feature = Feature::Padding(position, 13);
    ///
    /// assert_eq!(feature.position(), position);
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn position(&self) -> Position {
        match self {
            Self::Bases(pos, _) => *pos,
            Self::Scores(pos, _) => *pos,
            Self::ReadBase(pos, _, _) => *pos,
            Self::Substitution(pos, _) => *pos,
            Self::Insertion(pos, _) => *pos,
            Self::Deletion(pos, _) => *pos,
            Self::InsertBase(pos, _) => *pos,
            Self::QualityScore(pos, _) => *pos,
            Self::ReferenceSkip(pos, _) => *pos,
            Self::SoftClip(pos, _) => *pos,
            Self::Padding(pos, _) => *pos,
            Self::HardClip(pos, _) => *pos,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_code() {
        let position = Position::MIN;

        assert_eq!(Feature::Bases(position, Vec::new()).code(), Code::Bases);
        assert_eq!(Feature::Scores(position, Vec::new()).code(), Code::Scores);
        assert_eq!(
            Feature::ReadBase(position, Base::N, 0).code(),
            Code::ReadBase
        );
        assert_eq!(
            Feature::Substitution(position, substitution::Value::Code(0)).code(),
            Code::Substitution
        );
        assert_eq!(
            Feature::Insertion(position, Vec::new()).code(),
            Code::Insertion
        );
        assert_eq!(Feature::Deletion(position, 0).code(), Code::Deletion);
        assert_eq!(
            Feature::InsertBase(position, Base::N).code(),
            Code::InsertBase
        );
        assert_eq!(
            Feature::QualityScore(position, 0).code(),
            Code::QualityScore
        );
        assert_eq!(
            Feature::ReferenceSkip(position, 0).code(),
            Code::ReferenceSkip
        );
        assert_eq!(
            Feature::SoftClip(position, Vec::new()).code(),
            Code::SoftClip
        );
        assert_eq!(Feature::Padding(position, 0).code(), Code::Padding);
        assert_eq!(Feature::HardClip(position, 0).code(), Code::HardClip);
    }

    #[test]
    fn test_position() {
        let position = Position::MIN;

        assert_eq!(Feature::Bases(position, Vec::new()).position(), position);
        assert_eq!(Feature::Scores(position, Vec::new()).position(), position);
        assert_eq!(Feature::ReadBase(position, Base::N, 0).position(), position);
        assert_eq!(
            Feature::Substitution(position, substitution::Value::Code(0)).position(),
            position
        );
        assert_eq!(
            Feature::Insertion(position, Vec::new()).position(),
            position
        );
        assert_eq!(Feature::Deletion(position, 0).position(), position);
        assert_eq!(Feature::InsertBase(position, Base::N).position(), position);
        assert_eq!(Feature::QualityScore(position, 0).position(), position);
        assert_eq!(Feature::ReferenceSkip(position, 0).position(), position);
        assert_eq!(Feature::SoftClip(position, Vec::new()).position(), position);
        assert_eq!(Feature::Padding(position, 0).position(), position);
        assert_eq!(Feature::HardClip(position, 0).position(), position);
    }
}
