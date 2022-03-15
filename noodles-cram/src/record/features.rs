mod with_positions;

pub use self::with_positions::WithPositions;

use std::{
    ops::{Deref, DerefMut},
    slice,
};

use noodles_core::Position;
use noodles_sam as sam;

use super::Feature;

/// CRAM record features.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Features(Vec<Feature>);

impl Features {
    /// Converts SAM record CIGAR operations to CRAM record features.
    pub fn from_cigar(cigar: &sam::record::Cigar, sequence: &sam::record::Sequence) -> Self {
        cigar_to_features(cigar, sequence)
    }

    pub(crate) fn with_positions(
        &self,
        alignment_start: Position,
    ) -> WithPositions<'_, slice::Iter<'_, Feature>> {
        WithPositions::new(self.iter(), alignment_start)
    }
}

impl Deref for Features {
    type Target = Vec<Feature>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Features {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<Feature>> for Features {
    fn from(features: Vec<Feature>) -> Self {
        Self(features)
    }
}

fn cigar_to_features(cigar: &sam::record::Cigar, sequence: &sam::record::Sequence) -> Features {
    use sam::record::{cigar::op::Kind, quality_scores::Score};

    let mut features = Features::default();

    // SAFETY: 1 is non-zero.
    let mut read_position = Position::new(1).unwrap();

    for op in cigar.iter() {
        let feature = match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if op.len() == 1 {
                    let base = sequence[read_position];
                    // FIXME: quality score
                    let score = Score::default();
                    Feature::ReadBase(read_position, base, score)
                } else {
                    let end = read_position
                        .checked_add(op.len())
                        .expect("attempt to add with overflow");
                    let bases = sequence[read_position..end].to_vec();
                    Feature::Bases(read_position, bases)
                }
            }
            Kind::Insertion => {
                if op.len() == 1 {
                    Feature::InsertBase(read_position, sequence[read_position])
                } else {
                    let end = read_position
                        .checked_add(op.len())
                        .expect("attempt to add with overflow");
                    let bases = sequence[read_position..end].to_vec();
                    Feature::Insertion(read_position, bases)
                }
            }
            Kind::Deletion => Feature::Deletion(read_position, op.len()),
            Kind::Skip => Feature::ReferenceSkip(read_position, op.len()),
            Kind::SoftClip => {
                let end = read_position
                    .checked_add(op.len())
                    .expect("attempt to add with overflow");
                let bases = sequence[read_position..end].to_vec();
                Feature::SoftClip(read_position, bases)
            }
            Kind::HardClip => Feature::HardClip(read_position, op.len()),
            Kind::Pad => Feature::Padding(read_position, op.len()),
        };

        features.push(feature);

        if matches!(
            op.kind(),
            Kind::Match
                | Kind::Insertion
                | Kind::SoftClip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
        ) {
            read_position = read_position
                .checked_add(op.len())
                .expect("attempt to add with overflow");
        }
    }

    features
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_to_features() -> Result<(), Box<dyn std::error::Error>> {
        use sam::record::{quality_scores::Score, sequence::Base};

        let cigar = "1M".parse()?;
        let sequence = "A".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![Feature::ReadBase(
            Position::try_from(1)?,
            Base::A,
            Score::default(),
        )]);
        assert_eq!(actual, expected);

        let cigar = "2M".parse()?;
        let sequence = "AC".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![Feature::Bases(
            Position::try_from(1)?,
            vec![Base::A, Base::C],
        )]);
        assert_eq!(actual, expected);

        let cigar = "1I2M".parse()?;
        let sequence = "ACG".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![
            Feature::InsertBase(Position::try_from(1)?, Base::A),
            Feature::Bases(Position::try_from(2)?, vec![Base::C, Base::G]),
        ]);
        assert_eq!(actual, expected);

        let cigar = "2I2M".parse()?;
        let sequence = "ACGT".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![
            Feature::Insertion(Position::try_from(1)?, vec![Base::A, Base::C]),
            Feature::Bases(Position::try_from(3)?, vec![Base::G, Base::T]),
        ]);
        assert_eq!(actual, expected);

        let cigar = "1D2M".parse()?;
        let sequence = "AC".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![
            Feature::Deletion(Position::try_from(1)?, 1),
            Feature::Bases(Position::try_from(1)?, vec![Base::A, Base::C]),
        ]);
        assert_eq!(actual, expected);

        let cigar = "1N2M".parse()?;
        let sequence = "AC".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![
            Feature::ReferenceSkip(Position::try_from(1)?, 1),
            Feature::Bases(Position::try_from(1)?, vec![Base::A, Base::C]),
        ]);
        assert_eq!(actual, expected);

        let cigar = "1S2M".parse()?;
        let sequence = "ACG".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![
            Feature::SoftClip(Position::try_from(1)?, vec![Base::A]),
            Feature::Bases(Position::try_from(2)?, vec![Base::C, Base::G]),
        ]);
        assert_eq!(actual, expected);

        let cigar = "1H2M".parse()?;
        let sequence = "AC".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![
            Feature::HardClip(Position::try_from(1)?, 1),
            Feature::Bases(Position::try_from(1)?, vec![Base::A, Base::C]),
        ]);
        assert_eq!(actual, expected);

        let cigar = "1P2M".parse()?;
        let sequence = "AC".parse()?;
        let actual = cigar_to_features(&cigar, &sequence);
        let expected = Features::from(vec![
            Feature::Padding(Position::try_from(1)?, 1),
            Feature::Bases(Position::try_from(1)?, vec![Base::A, Base::C]),
        ]);
        assert_eq!(actual, expected);

        Ok(())
    }
}
