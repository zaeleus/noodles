use std::io;

use bstr::BStr;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::data::field::{Tag, Value},
        record_buf::{data::field::Value as ValueBuf, QualityScores, Sequence},
    },
};

use super::{Feature, Record};
use crate::record::{Flags, MateFlags};

impl Record {
    pub fn try_from_alignment_record(
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<Self> {
        let bam_flags = record.flags()?;
        let mut cram_flags = Flags::default();

        let sequence = Sequence::from(record.sequence().iter().collect::<Vec<_>>());

        let quality_scores = if record.quality_scores().is_empty() {
            QualityScores::default()
        } else {
            if bam_flags.is_unmapped() {
                cram_flags.insert(Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY);
            }

            QualityScores::from(
                record
                    .quality_scores()
                    .iter()
                    .collect::<io::Result<Vec<_>>>()?,
            )
        };

        let features = cigar_to_features(
            record.cigar().as_ref(),
            cram_flags,
            &sequence,
            &quality_scores,
        )?;

        let data = record.data();
        let (data_buf, read_group_name) = get_filtered_data(data.as_ref())?;

        let read_group_id = read_group_name
            .map(|name| get_read_group_id(header, name))
            .transpose()?;

        Ok(Self {
            bam_flags,
            cram_flags,
            reference_sequence_id: record.reference_sequence_id(header).transpose()?,
            read_length: record.sequence().len(),
            alignment_start: record.alignment_start().transpose()?,
            read_group_id,
            name: record.name().map(|s| s.into()),
            mate_flags: MateFlags::default(),
            mate_reference_sequence_id: record.mate_reference_sequence_id(header).transpose()?,
            mate_alignment_start: record.mate_alignment_start().transpose()?,
            template_length: record.template_length()?,
            mate_distance: None,
            data: data_buf,
            features,
            mapping_quality: record.mapping_quality().transpose()?,
            sequence: sequence.into(),
            quality_scores: quality_scores.into(),
        })
    }
}

fn get_read_group_id(header: &sam::Header, read_group_name: &BStr) -> io::Result<usize> {
    header
        .read_groups()
        .get_index_of(read_group_name)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid read group name"))
}

fn cigar_to_features(
    cigar: &dyn sam::alignment::record::Cigar,
    flags: Flags,
    sequence: &Sequence,
    quality_scores: &QualityScores,
) -> io::Result<Vec<Feature>> {
    use noodles_sam::alignment::record::cigar::op::Kind;

    let mut features = Vec::new();
    let mut position = Position::MIN;

    for result in cigar.iter() {
        let op = result?;

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

                    if !flags.quality_scores_are_stored_as_array() {
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

                    if !flags.quality_scores_are_stored_as_array() {
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

                    if !flags.quality_scores_are_stored_as_array() {
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

                if !flags.quality_scores_are_stored_as_array() {
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
        }

        if op.kind().consumes_read() {
            position = position
                .checked_add(op.len())
                .expect("attempt to add with overflow");
        }
    }

    Ok(features)
}

#[allow(clippy::type_complexity)]
fn get_filtered_data(
    data: &dyn sam::alignment::record::Data,
) -> io::Result<(Vec<(Tag, ValueBuf)>, Option<&BStr>)> {
    let mut data_buf = Vec::new();
    let mut read_group_name = None;

    for result in data.iter() {
        let (tag, value) = result?;

        if tag == Tag::READ_GROUP {
            let Value::String(s) = value else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid read group field value",
                ));
            };

            read_group_name = Some(s);
        }

        data_buf.push((tag, value.try_into()?));
    }

    Ok((data_buf, read_group_name))
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::{
        record::cigar::{op::Kind, Op},
        record_buf::Cigar,
    };

    use super::*;

    #[test]
    fn test_cigar_to_features() -> Result<(), Box<dyn std::error::Error>> {
        let flags = Flags::default();

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = QualityScores::from(vec![45]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![Feature::ReadBase {
            position: Position::try_from(1)?,
            base: b'A',
            quality_score: 45,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = QualityScores::from(vec![45, 35]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::Bases {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = QualityScores::from(vec![45, 35]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::InsertBase {
                position: Position::try_from(1)?,
                base: b'A',
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::ReadBase {
                position: Position::try_from(2)?,
                base: b'C',
                quality_score: 35,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = QualityScores::from(vec![45, 35, 43]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::ReadBase {
                position: Position::try_from(3)?,
                base: b'G',
                quality_score: 43,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = QualityScores::from(vec![45, 35]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::Deletion {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Bases {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Skip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = QualityScores::from(vec![45]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::ReferenceSkip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::ReadBase {
                position: Position::try_from(1)?,
                base: b'A',
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = QualityScores::from(vec![45, 35]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::ReadBase {
                position: Position::try_from(2)?,
                base: b'C',
                quality_score: 35,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = QualityScores::from(vec![45, 35, 43]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::ReadBase {
                position: Position::try_from(3)?,
                base: b'G',
                quality_score: 43,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = QualityScores::from(vec![45]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::HardClip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::ReadBase {
                position: Position::try_from(1)?,
                base: b'A',
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Pad, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = QualityScores::from(vec![45]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::Padding {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::ReadBase {
                position: Position::try_from(1)?,
                base: b'A',
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_cigar_to_features_with_quality_scores_stored_as_array(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let flags = Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = QualityScores::from(vec![45]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![Feature::ReadBase {
            position: Position::try_from(1)?,
            base: b'A',
            quality_score: 45,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = QualityScores::from(vec![45, 35]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![Feature::Bases {
            position: Position::try_from(1)?,
            bases: vec![b'A', b'C'],
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = QualityScores::from(vec![45, 35]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::InsertBase {
                position: Position::try_from(1)?,
                base: b'A',
            },
            Feature::ReadBase {
                position: Position::try_from(2)?,
                base: b'C',
                quality_score: 35,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = QualityScores::from(vec![45, 35, 43]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::ReadBase {
                position: Position::try_from(3)?,
                base: b'G',
                quality_score: 43,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = QualityScores::from(vec![45, 35]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::Deletion {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Bases {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Skip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = QualityScores::from(vec![45]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::ReferenceSkip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::ReadBase {
                position: Position::try_from(1)?,
                base: b'A',
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = QualityScores::from(vec![45, 35]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::ReadBase {
                position: Position::try_from(2)?,
                base: b'C',
                quality_score: 35,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = QualityScores::from(vec![45, 35, 43]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::ReadBase {
                position: Position::try_from(3)?,
                base: b'G',
                quality_score: 43,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = QualityScores::from(vec![45]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::HardClip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::ReadBase {
                position: Position::try_from(1)?,
                base: b'A',
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Pad, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = QualityScores::from(vec![45]);
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![
            Feature::Padding {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::ReadBase {
                position: Position::try_from(1)?,
                base: b'A',
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        Ok(())
    }
}
