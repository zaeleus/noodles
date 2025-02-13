use std::{io, slice};

use noodles_core::Position;

use crate::record::Feature;

const MISSING: u8 = 0;

pub(super) struct Iter<'r> {
    features: slice::Iter<'r, Feature>,
    read_length: usize,
    read_position: Position,
    state: State<'r>,
}

impl<'r> Iter<'r> {
    pub(super) fn new(features: &'r [Feature], read_length: usize) -> Self {
        Self {
            features: features.iter(),
            read_length,
            read_position: Position::MIN,
            state: State::Next,
        }
    }
}

enum State<'r> {
    Next,
    Missing {
        end: Position,
        next_feature: &'r Feature,
    },
    Prepare(&'r Feature),
    Score(u8),
    Scores(slice::Iter<'r, u8>),
    Finish,
    Done,
}

impl Iterator for Iter<'_> {
    type Item = io::Result<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        fn succ(position: Position) -> Position {
            position
                .checked_add(1)
                .expect("attempt to add with overflow")
        }

        loop {
            match self.state {
                State::Next => {
                    self.state = if let Some(feature) = self.features.next() {
                        let position = feature.position();

                        if self.read_position < position {
                            State::Missing {
                                end: position,
                                next_feature: feature,
                            }
                        } else {
                            State::Prepare(feature)
                        }
                    } else if usize::from(self.read_position) < self.read_length {
                        State::Finish
                    } else {
                        State::Done
                    };
                }
                State::Missing { end, next_feature } => {
                    if self.read_position < end {
                        self.read_position = succ(self.read_position);
                        return Some(Ok(MISSING));
                    } else {
                        self.state = State::Prepare(next_feature);
                    }
                }
                State::Prepare(feature) => {
                    self.state = match feature {
                        Feature::Scores { quality_scores, .. } => {
                            State::Scores(quality_scores.iter())
                        }
                        Feature::ReadBase { quality_score, .. }
                        | Feature::QualityScore { quality_score, .. } => {
                            State::Score(*quality_score)
                        }
                        _ => State::Next,
                    };
                }
                State::Score(score) => {
                    self.read_position = succ(self.read_position);
                    self.state = State::Next;
                    return Some(Ok(score));
                }
                State::Scores(ref mut iter) => match iter.next() {
                    Some(score) => {
                        self.read_position = succ(self.read_position);
                        return Some(Ok(*score));
                    }
                    None => self.state = State::Next,
                },
                State::Finish => {
                    if usize::from(self.read_position) <= self.read_length {
                        self.read_position = succ(self.read_position);
                        return Some(Ok(MISSING));
                    } else {
                        self.state = State::Done;
                    }
                }
                State::Done => return None,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() -> Result<(), Box<dyn std::error::Error>> {
        let features = [
            Feature::ReadBase {
                position: Position::try_from(2)?,
                base: b'A',
                quality_score: 5,
            },
            Feature::QualityScore {
                position: Position::try_from(3)?,
                quality_score: 8,
            },
            Feature::Scores {
                position: Position::try_from(5)?,
                quality_scores: vec![13, 21],
            },
        ];

        let actual: Vec<_> = Iter::new(&features, 8).collect::<io::Result<_>>()?;
        assert_eq!(actual, [MISSING, 5, 8, MISSING, 13, 21, MISSING, MISSING]);

        Ok(())
    }
}
