mod cigar;
mod with_positions;

pub use self::{cigar::Cigar, with_positions::WithPositions};

use std::{
    io,
    ops::{Deref, DerefMut},
    slice,
};

use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::record_buf::{QualityScores, Sequence},
};

use super::{Feature, Flags};

/// CRAM record features.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Features(Vec<Feature>);

impl Features {
    /// Converts SAM record CIGAR operations to CRAM record features.
    pub fn from_cigar(
        flags: Flags,
        cigar: &sam::alignment::record_buf::Cigar,
        sequence: &Sequence,
        quality_scores: &QualityScores,
    ) -> Self {
        cigar_to_features(flags, cigar, sequence, quality_scores)
    }

    /// Converts CRAM features to SAM CIGAR operations.
    pub fn try_into_cigar(
        &self,
        read_length: usize,
    ) -> io::Result<sam::alignment::record_buf::Cigar> {
        use sam::alignment::record::cigar::{op::Kind, Op};

        fn merge_or_insert_op(ops: &mut Vec<(Kind, usize)>, kind: Kind, len: usize) {
            if let Some(last_op) = ops.last_mut() {
                if last_op.0 == kind {
                    last_op.1 += len;
                    return;
                }
            }

            ops.push((kind, len));
        }

        let mut ops = Vec::new();
        let mut read_position = Position::MIN;

        for feature in self.iter() {
            if feature.position() > read_position {
                let len = usize::from(feature.position()) - usize::from(read_position);
                merge_or_insert_op(&mut ops, Kind::Match, len);
                read_position = feature.position();
            }

            let (kind, len) = match feature {
                Feature::Substitution { .. } => (Kind::Match, 1),
                Feature::Insertion { bases, .. } => (Kind::Insertion, bases.len()),
                Feature::Deletion { len, .. } => (Kind::Deletion, *len),
                Feature::InsertBase { .. } => (Kind::Insertion, 1),
                Feature::ReferenceSkip { len, .. } => (Kind::Skip, *len),
                Feature::SoftClip { bases, .. } => (Kind::SoftClip, bases.len()),
                Feature::Padding { len, .. } => (Kind::Pad, *len),
                Feature::HardClip { len, .. } => (Kind::HardClip, *len),
                _ => continue,
            };

            merge_or_insert_op(&mut ops, kind, len);

            if kind.consumes_read() {
                read_position = read_position
                    .checked_add(len)
                    .expect("attempt to add with overflow");
            }
        }

        if usize::from(read_position) <= read_length {
            let len = read_length - usize::from(read_position) + 1;
            merge_or_insert_op(&mut ops, Kind::Match, len);
        }

        Ok(ops
            .into_iter()
            .map(|(kind, len)| Op::new(kind, len))
            .collect())
    }

    /// Returns an iterator over features as CIGAR operations.
    pub fn cigar(&self, read_length: usize) -> Cigar<'_> {
        Cigar::new(&self.0, read_length)
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

fn cigar_to_features(
    flags: Flags,
    cigar: &sam::alignment::record_buf::Cigar,
    sequence: &Sequence,
    quality_scores: &QualityScores,
) -> Features {
    use sam::alignment::record::cigar::op::Kind;

    let mut features = Features::default();
    let mut position = Position::MIN;

    for op in cigar.as_ref().iter() {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if op.len() == 1 {
                    let base = sequence[position];
                    let quality_score = quality_scores[position];

                    features.push(Feature::ReadBase {
                        position,
                        base,
                        quality_score,
                    });
                } else {
                    let end = position
                        .checked_add(op.len())
                        .expect("attempt to add with overflow");

                    let bases = sequence[position..end].to_vec();
                    features.push(Feature::Bases { position, bases });

                    if !flags.are_quality_scores_stored_as_array() {
                        let quality_scores = quality_scores[position..end].to_vec();

                        features.push(Feature::Scores {
                            position,
                            quality_scores,
                        });
                    }
                }
            }
            Kind::Insertion => {
                if op.len() == 1 {
                    let base = sequence[position];
                    features.push(Feature::InsertBase { position, base });

                    if !flags.are_quality_scores_stored_as_array() {
                        let quality_score = quality_scores[position];

                        features.push(Feature::QualityScore {
                            position,
                            quality_score,
                        });
                    }
                } else {
                    let end = position
                        .checked_add(op.len())
                        .expect("attempt to add with overflow");

                    let bases = sequence[position..end].to_vec();
                    features.push(Feature::Insertion { position, bases });

                    if !flags.are_quality_scores_stored_as_array() {
                        let quality_scores = quality_scores[position..end].to_vec();

                        features.push(Feature::Scores {
                            position,
                            quality_scores,
                        });
                    }
                }
            }
            Kind::Deletion => features.push(Feature::Deletion {
                position,
                len: op.len(),
            }),
            Kind::Skip => features.push(Feature::ReferenceSkip {
                position,
                len: op.len(),
            }),
            Kind::SoftClip => {
                let end = position
                    .checked_add(op.len())
                    .expect("attempt to add with overflow");

                let bases = &sequence[position..end];

                features.push(Feature::SoftClip {
                    position,
                    bases: bases.to_vec(),
                });

                if !flags.are_quality_scores_stored_as_array() {
                    if bases.len() == 1 {
                        let quality_score = quality_scores[position];

                        features.push(Feature::QualityScore {
                            position,
                            quality_score,
                        });
                    } else {
                        let quality_scores = quality_scores[position..end].to_vec();

                        features.push(Feature::Scores {
                            position,
                            quality_scores,
                        });
                    }
                }
            }
            Kind::HardClip => features.push(Feature::HardClip {
                position,
                len: op.len(),
            }),
            Kind::Pad => features.push(Feature::Padding {
                position,
                len: op.len(),
            }),
        };

        if op.kind().consumes_read() {
            position = position
                .checked_add(op.len())
                .expect("attempt to add with overflow");
        }
    }

    features
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::cigar::{op::Kind, Op};

    use super::*;

    #[test]
    fn test_try_into_cigar() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::feature::substitution;

        let features = Features::default();
        assert_eq!(
            features.try_into_cigar(4)?,
            [Op::new(Kind::Match, 4)].into_iter().collect()
        );

        let features = Features::from(vec![Feature::SoftClip {
            position: Position::try_from(1)?,
            bases: vec![b'A', b'T'],
        }]);
        assert_eq!(
            features.try_into_cigar(4)?,
            [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 2)]
                .into_iter()
                .collect()
        );

        let features = Features::from(vec![Feature::SoftClip {
            position: Position::try_from(4)?,
            bases: vec![b'G'],
        }]);
        assert_eq!(
            features.try_into_cigar(4)?,
            [Op::new(Kind::Match, 3), Op::new(Kind::SoftClip, 1)]
                .into_iter()
                .collect()
        );

        let features = Features::from(vec![Feature::HardClip {
            position: Position::try_from(1)?,
            len: 2,
        }]);
        assert_eq!(
            features.try_into_cigar(4)?,
            [Op::new(Kind::HardClip, 2), Op::new(Kind::Match, 4)]
                .into_iter()
                .collect()
        );

        let features = Features::from(vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                value: substitution::Value::Code(0),
            },
        ]);
        assert_eq!(
            features.try_into_cigar(4)?,
            [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 3)]
                .into_iter()
                .collect()
        );

        let features = Features::from(vec![Feature::Substitution {
            position: Position::try_from(2)?,
            value: substitution::Value::Code(0),
        }]);
        assert_eq!(
            features.try_into_cigar(4)?,
            [Op::new(Kind::Match, 4)].into_iter().collect()
        );

        Ok(())
    }
}
