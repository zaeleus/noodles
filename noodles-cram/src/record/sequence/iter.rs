mod with_positions;

use std::slice;

use noodles_core::Position;
use noodles_fasta as fasta;

use self::with_positions::WithPositions;
use crate::{
    container::compression_header::preservation_map::{
        substitution_matrix::Base, SubstitutionMatrix,
    },
    record::Feature,
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
                                let raw_reference_base = reference_sequence[reference_position];

                                let reference_base = Base::try_from(raw_reference_base)
                                    .expect("invalid reference base");

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
        (self.read_length, Some(self.read_length))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"CGTCCGTAACACTAGG".to_vec());

        let features = [
            Feature::Bases {
                position: Position::try_from(2)?,
                bases: b"A",
            },
            Feature::Scores {
                position: Position::try_from(2)?,
                quality_scores: &[0],
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
            Feature::Insertion {
                position: Position::try_from(5)?,
                bases: b"G",
            },
            Feature::Deletion {
                position: Position::try_from(6)?,
                len: 1,
            },
            Feature::InsertBase {
                position: Position::try_from(6)?,
                base: b'A',
            },
            Feature::QualityScore {
                position: Position::try_from(6)?,
                quality_score: 0,
            },
            Feature::ReferenceSkip {
                position: Position::try_from(7)?,
                len: 1,
            },
            Feature::SoftClip {
                position: Position::try_from(7)?,
                bases: b"T",
            },
            Feature::Padding {
                position: Position::try_from(8)?,
                len: 1,
            },
            Feature::HardClip {
                position: Position::try_from(8)?,
                len: 1,
            },
        ];

        let iter = Iter::new(
            Some(&reference_sequence),
            SubstitutionMatrix::default(),
            &features,
            Position::MIN,
            8,
        );

        let actual: Vec<_> = iter.collect();
        assert_eq!(actual, b"CACAGATT");

        Ok(())
    }

    #[test]
    fn test_size_hint() {
        let iter = Iter::new(None, SubstitutionMatrix::default(), &[], Position::MIN, 4);
        assert_eq!(iter.size_hint(), (4, Some(4)));
    }
}
