use std::slice;

use noodles_core::Position;

use crate::record::Feature;

pub(super) struct WithPositions<'r, 'c: 'r> {
    features: slice::Iter<'r, Feature<'c>>,
    reference_position: Position,
    read_position: Position,
}

impl<'r, 'c: 'r> WithPositions<'r, 'c> {
    pub(super) fn new(features: &'r [Feature<'c>], alignment_start: Position) -> Self {
        Self {
            features: features.iter(),
            reference_position: alignment_start,
            read_position: Position::MIN,
        }
    }

    pub(super) fn positions(&self) -> (Position, Position) {
        (self.reference_position, self.read_position)
    }
}

impl<'r, 'c: 'r> Iterator for WithPositions<'r, 'c> {
    type Item = ((Position, Position), &'r Feature<'c>);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let feature = self.features.next()?;

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
    use std::borrow::Cow;

    use super::*;

    #[test]
    fn test_next() {
        let features = [
            Feature::Bases {
                position: Position::MIN,
                bases: Cow::Borrowed(b"AC"),
            },
            Feature::Scores {
                position: Position::MIN,
                quality_scores: Cow::Borrowed(&[0, 0]),
            },
        ];

        let mut iter = WithPositions::new(&features, Position::MIN);

        assert_eq!(
            iter.next(),
            Some(((Position::MIN, Position::MIN), &features[0]))
        );
        assert!(iter.next().is_none());
    }
}
