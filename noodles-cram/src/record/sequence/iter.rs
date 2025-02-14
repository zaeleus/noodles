use std::slice;

use noodles_core::Position;
use noodles_fasta as fasta;

use crate::{
    container::compression_header::preservation_map::SubstitutionMatrix,
    record::{feature::substitution::Base, Feature},
};

pub(super) struct Iter<'r, 'c: 'r> {
    reference_sequence: Option<&'c fasta::record::Sequence>,
    substitution_matrix: SubstitutionMatrix,
    features: WithPositions<'r, 'c>,
    read_length: usize,
    last_reference_position: Position,
    last_read_position: Position,
    state: State<'c>,
}

impl<'r, 'c: 'r> Iter<'r, 'c> {
    pub(super) fn new(
        reference_sequence: Option<&'c fasta::record::Sequence>,
        substitution_matrix: SubstitutionMatrix,
        features: &'r [Feature<'c>],
        alignment_start: Position,
        read_length: usize,
    ) -> Self {
        let features = WithPositions::new(features, alignment_start);
        let (last_reference_position, last_read_position) = features.positions();

        Self {
            reference_sequence,
            substitution_matrix,
            features,
            read_length,
            last_reference_position,
            last_read_position,
            state: State::Next,
        }
    }
}

enum State<'c> {
    Next,
    Prepare(slice::Iter<'c, u8>, Position, &'c Feature<'c>),
    Base(u8),
    Bases(slice::Iter<'c, u8>),
    Finish(slice::Iter<'c, u8>),
    Done,
}

impl<'r: 'c, 'c: 'r> Iterator for Iter<'r, 'c> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state {
                State::Next => {
                    self.state = if let Some(((reference_position, read_position), feature)) =
                        self.features.next()
                    {
                        let bases = if let Some(reference_sequence) = self.reference_sequence {
                            &reference_sequence[self.last_reference_position..reference_position]
                        } else if read_position != self.last_read_position {
                            panic!("next: missing reference sequence");
                        } else {
                            &[]
                        };

                        State::Prepare(bases.iter(), reference_position, feature)
                    } else if let Some(reference_sequence) = self.reference_sequence {
                        let last_read_position = usize::from(self.last_read_position);

                        if last_read_position > self.read_length {
                            State::Done
                        } else {
                            let len = self.read_length - last_read_position + 1;

                            let end = self
                                .last_reference_position
                                .checked_add(len)
                                .expect("attempt to add with overflow");

                            let bases = &reference_sequence[self.last_reference_position..end];

                            State::Finish(bases.iter())
                        }
                    } else if usize::from(self.last_read_position) != self.read_length + 1 {
                        panic!(
                            "next: missing reference sequence: {} != {}",
                            usize::from(self.last_read_position),
                            self.read_length
                        );
                    } else {
                        State::Done
                    };

                    let (next_reference_position, next_read_position) = self.features.positions();
                    self.last_reference_position = next_reference_position;
                    self.last_read_position = next_read_position;
                }
                State::Prepare(ref mut reference_bases, reference_position, feature) => {
                    if let Some(base) = reference_bases.next() {
                        return Some(*base);
                    }

                    self.state = match feature {
                        Feature::Bases { bases, .. } => State::Bases(bases.iter()),
                        Feature::Scores { .. } => State::Next,
                        Feature::ReadBase { base, .. } => State::Base(*base),
                        Feature::Substitution { code, .. } => {
                            if let Some(reference_sequence) = self.reference_sequence {
                                let reference_base =
                                    Base::try_from(reference_sequence[reference_position])
                                        .unwrap_or_default();

                                let read_base = self.substitution_matrix.get(reference_base, *code);

                                State::Base(u8::from(read_base))
                            } else {
                                panic!("missing reference sequence (substitution)");
                            }
                        }
                        Feature::Insertion { bases, .. } => State::Bases(bases.iter()),
                        Feature::Deletion { .. } => State::Next,
                        Feature::InsertBase { base, .. } => State::Base(*base),
                        Feature::QualityScore { .. } => State::Next,
                        Feature::ReferenceSkip { .. } => State::Next,
                        Feature::SoftClip { bases, .. } => State::Bases(bases.iter()),
                        Feature::Padding { .. } => State::Next,
                        Feature::HardClip { .. } => State::Next,
                    }
                }
                State::Base(base) => {
                    self.state = State::Next;
                    return Some(base);
                }
                State::Bases(ref mut iter) => match iter.next() {
                    Some(base) => return Some(*base),
                    None => self.state = State::Next,
                },
                State::Finish(ref mut iter) => match iter.next() {
                    Some(base) => return Some(*base),
                    None => self.state = State::Done,
                },
                State::Done => return None,
            }
        }
    }
}

struct WithPositions<'r, 'c: 'r> {
    features: slice::Iter<'r, Feature<'c>>,
    reference_position: Position,
    read_position: Position,
}

impl<'r, 'c: 'r> WithPositions<'r, 'c> {
    fn new(features: &'r [Feature<'c>], alignment_start: Position) -> Self {
        Self {
            features: features.iter(),
            reference_position: alignment_start,
            read_position: Position::MIN,
        }
    }

    fn positions(&self) -> (Position, Position) {
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

        let mut iter = WithPositions::new(&features, Position::MIN);

        assert_eq!(
            iter.next(),
            Some(((Position::MIN, Position::MIN), &features[0]))
        );
        assert!(iter.next().is_none());

        Ok(())
    }
}
