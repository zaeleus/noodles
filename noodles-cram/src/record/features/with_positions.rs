use noodles_core::Position;

use crate::record::Feature;

pub struct WithPositions<'r, 'c: 'r, I>
where
    I: Iterator<Item = &'r Feature<'c>>,
{
    iter: I,
    reference_position: Position,
    read_position: Position,
}

impl<'r, 'c: 'r, I> WithPositions<'r, 'c, I>
where
    I: Iterator<Item = &'r Feature<'c>>,
{
    pub fn new(iter: I, alignment_start: Position) -> Self {
        Self {
            iter,
            reference_position: alignment_start,
            read_position: Position::MIN,
        }
    }

    /// Returns the current reference position and read position.
    ///
    /// These are 1-based.
    pub fn positions(&self) -> (Position, Position) {
        (self.reference_position, self.read_position)
    }
}

impl<'r, 'c: 'r, I> Iterator for WithPositions<'r, 'c, I>
where
    I: Iterator<Item = &'r Feature<'c>>,
{
    type Item = ((Position, Position), I::Item);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let feature = self.iter.next()?;

            let (reference_position_delta, read_position_delta) = match feature {
                Feature::Bases { bases, .. } => (bases.len(), bases.len()),
                Feature::Scores { .. } => continue,
                Feature::ReadBase { .. } => (1, 1),
                Feature::Substitution { .. } => (1, 1),
                Feature::Insertion { bases, .. } => (0, bases.len()),
                Feature::Deletion { len, .. } => (*len, 0),
                Feature::InsertBase { .. } => (0, 1),
                Feature::QualityScore { .. } => continue,
                Feature::ReferenceSkip { len, .. } => (*len, 0),
                Feature::SoftClip { bases, .. } => (0, bases.len()),
                Feature::Padding { .. } => (0, 0),
                Feature::HardClip { .. } => (0, 0),
            };

            let feature_position = usize::from(feature.position());
            let match_len = feature_position - usize::from(self.read_position);

            self.reference_position = self
                .reference_position
                .checked_add(match_len)
                .expect("attempt to add with overflow");

            self.read_position = self
                .read_position
                .checked_add(match_len)
                .expect("attempt to add with overflow");

            let positions = self.positions();

            self.reference_position = self
                .reference_position
                .checked_add(reference_position_delta)
                .expect("attempt to add with overflow");

            self.read_position = self
                .read_position
                .checked_add(read_position_delta)
                .expect("attempt to add with overflow");

            return Some((positions, feature));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let features = [
            Feature::Bases {
                position: Position::MIN,
                bases: b"AC",
            },
            Feature::Scores {
                position: Position::MIN,
                quality_scores: &[0, 0],
            },
        ];

        let mut iter = WithPositions::new(features.iter(), Position::MIN);

        assert_eq!(
            iter.next(),
            Some(((Position::MIN, Position::MIN), &features[0]))
        );
        assert!(iter.next().is_none());

        Ok(())
    }
}
