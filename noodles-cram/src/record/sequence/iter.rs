mod with_positions;

use std::{iter::FusedIterator, ops::Range, slice};

use noodles_core::Position;

use self::with_positions::WithPositions;
use crate::{
    container::compression_header::preservation_map::{
        SubstitutionMatrix, substitution_matrix::Base,
    },
    record::Feature,
};

struct ReferenceSequence<'c>(&'c [u8]);

impl<'c> ReferenceSequence<'c> {
    fn new(src: &'c [u8]) -> Self {
        Self(src)
    }

    fn at(&self, position: Position) -> u8 {
        let i = usize::from(position) - 1;
        self.0[i]
    }

    fn slice(&self, range: Range<Position>) -> &'c [u8] {
        let start = usize::from(range.start) - 1;
        let end = usize::from(range.end) - 1;
        &self.0[start..end]
    }
}

pub(super) struct Iter<'r, 'c: 'r> {
    reference_sequence: Option<ReferenceSequence<'c>>,
    substitution_matrix: SubstitutionMatrix,
    features: WithPositions<'r, 'c>,
    read_length: usize,
    last_reference_position: Position,
    last_read_position: Position,
    state: State<'r, 'c>,
}

impl<'r, 'c: 'r> Iter<'r, 'c> {
    pub(super) fn new(
        reference_sequence: Option<&'c [u8]>,
        substitution_matrix: SubstitutionMatrix,
        features: &'r [Feature<'c>],
        alignment_start: Position,
        read_length: usize,
    ) -> Self {
        let features = WithPositions::new(features, alignment_start);
        let (last_reference_position, last_read_position) = features.positions();

        Self {
            reference_sequence: reference_sequence.map(ReferenceSequence::new),
            substitution_matrix,
            features,
            read_length,
            last_reference_position,
            last_read_position,
            state: State::Next,
        }
    }
}

enum State<'r, 'c: 'r> {
    Next,
    Prepare(slice::Iter<'c, u8>, Position, &'r Feature<'c>),
    Base(u8),
    Bases(slice::Iter<'r, u8>),
    Finish(slice::Iter<'c, u8>),
    Done,
}

impl<'r, 'c: 'r> Iterator for Iter<'r, 'c> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state {
                State::Next => {
                    self.state = if let Some(((reference_position, read_position), feature)) =
                        self.features.next()
                    {
                        let bases =
                            if let Some(reference_sequence) = self.reference_sequence.as_ref() {
                                reference_sequence
                                    .slice(self.last_reference_position..reference_position)
                            } else if read_position != self.last_read_position {
                                panic!("next: missing reference sequence");
                            } else {
                                &[]
                            };

                        State::Prepare(bases.iter(), reference_position, feature)
                    } else if let Some(reference_sequence) = self.reference_sequence.as_ref() {
                        let last_read_position = usize::from(self.last_read_position);

                        if last_read_position > self.read_length {
                            State::Done
                        } else {
                            let len = self.read_length - last_read_position + 1;

                            let end = self
                                .last_reference_position
                                .checked_add(len)
                                .expect("attempt to add with overflow");

                            let bases = reference_sequence.slice(self.last_reference_position..end);

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
                            if let Some(reference_sequence) = self.reference_sequence.as_ref() {
                                let raw_reference_base = reference_sequence.at(reference_position);

                                let reference_base =
                                    Base::try_from(raw_reference_base).unwrap_or(Base::N);

                                let read_base = self.substitution_matrix.get(reference_base, *code);

                                let raw_read_base = if raw_reference_base.is_ascii_lowercase() {
                                    u8::from(read_base).to_ascii_lowercase()
                                } else {
                                    u8::from(read_base)
                                };

                                State::Base(raw_read_base)
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

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.read_length - (usize::from(self.last_read_position) - 1);

        match &self.state {
            State::Next => (n, Some(n)),
            State::Prepare(iter, ..) => {
                let m = n + iter.len() + 1;
                (m, Some(m))
            }
            State::Base(_) => {
                let m = n + 1;
                (m, Some(m))
            }
            State::Bases(iter) => {
                let m = n + iter.len();
                (m, Some(m))
            }
            State::Finish(iter) => iter.size_hint(),
            State::Done => (0, Some(0)),
        }
    }
}

impl<'r: 'c, 'c: 'r> ExactSizeIterator for Iter<'r, 'c> {}

impl<'r: 'c, 'c: 'r> FusedIterator for Iter<'r, 'c> {}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use super::*;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = b"CGTCYCGTAACACTAG";

        let features = [
            Feature::Bases {
                position: Position::try_from(2)?,
                bases: Cow::Borrowed(b"A"),
            },
            Feature::Scores {
                position: Position::try_from(2)?,
                quality_scores: Cow::Borrowed(&[0]),
            },
            Feature::ReadBase {
                position: Position::try_from(3)?,
                base: b'C',
                quality_score: 0,
            },
            Feature::Substitution {
                position: Position::try_from(4)?,
                code: 0b00,
            },
            Feature::Substitution {
                position: Position::try_from(5)?,
                code: 0b11,
            },
            Feature::Insertion {
                position: Position::try_from(6)?,
                bases: Cow::Borrowed(b"G"),
            },
            Feature::Deletion {
                position: Position::try_from(7)?,
                len: 1,
            },
            Feature::InsertBase {
                position: Position::try_from(7)?,
                base: b'A',
            },
            Feature::QualityScore {
                position: Position::try_from(7)?,
                quality_score: 0,
            },
            Feature::ReferenceSkip {
                position: Position::try_from(8)?,
                len: 1,
            },
            Feature::SoftClip {
                position: Position::try_from(8)?,
                bases: Cow::Borrowed(b"T"),
            },
            Feature::Padding {
                position: Position::try_from(9)?,
                len: 1,
            },
            Feature::HardClip {
                position: Position::try_from(9)?,
                len: 1,
            },
        ];

        let iter = Iter::new(
            Some(reference_sequence),
            SubstitutionMatrix::default(),
            &features,
            Position::MIN,
            9,
        );

        let actual: Vec<_> = iter.collect();
        assert_eq!(actual, b"CACATGATT");

        Ok(())
    }

    #[test]
    fn test_size_hint() {
        let iter = Iter::new(None, SubstitutionMatrix::default(), &[], Position::MIN, 4);
        assert_eq!(iter.size_hint(), (4, Some(4)));
    }
}
