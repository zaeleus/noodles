use std::io;

use bstr::BStr;
use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{
        record::data::field::{Tag, Value},
        record_buf::{QualityScores, Sequence, data::field::Value as ValueBuf},
    },
};

use super::{Feature, Record};
use crate::{
    container::compression_header::preservation_map::substitution_matrix::Base,
    record::{Flags, MateFlags},
};

impl Record {
    pub fn try_from_alignment_record(
        reference_sequence_repository: &fasta::Repository,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<Self> {
        let bam_flags = record.flags()?;
        let mut cram_flags = Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;

        let sequence = if record.sequence().is_empty() {
            cram_flags.insert(Flags::SEQUENCE_IS_MISSING);
            Sequence::default()
        } else {
            Sequence::from(record.sequence().iter().collect::<Vec<_>>())
        };

        let quality_scores = if record.quality_scores().is_empty() {
            QualityScores::default()
        } else {
            if bam_flags.is_unmapped() {
                cram_flags.insert(Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY);
            }

            record.quality_scores().iter().collect::<io::Result<_>>()?
        };

        let reference_sequence_id = record.reference_sequence_id(header).transpose()?;
        let alignment_start = record.alignment_start().transpose()?;

        let features = if let (Some(id), Some(start)) = (reference_sequence_id, alignment_start) {
            let (reference_sequence_name, _) = header
                .reference_sequences()
                .get_index(id)
                .expect("missing reference sequence ID");

            let reference_sequence = reference_sequence_repository
                .get(reference_sequence_name)
                .transpose()?
                .expect("missing reference sequence");

            cigar_to_features(
                record.cigar().as_ref(),
                reference_sequence.as_ref(),
                cram_flags,
                start,
                &sequence,
                &quality_scores,
            )
            .map(Some)?
        } else {
            None
        };

        let data = record.data();
        let (data_buf, read_group_name) = get_filtered_data(data.as_ref())?;

        let read_group_id = read_group_name
            .map(|name| get_read_group_id(header, name))
            .transpose()?;

        Ok(Self {
            bam_flags,
            cram_flags,
            reference_sequence_id,
            read_length: record.sequence().len(),
            alignment_start,
            read_group_id,
            name: record.name().map(|s| s.into()),
            mate_flags: MateFlags::default(),
            mate_reference_sequence_id: record.mate_reference_sequence_id(header).transpose()?,
            mate_alignment_start: record.mate_alignment_start().transpose()?,
            template_length: record.template_length()?,
            mate_distance: None,
            data: data_buf,
            features: features.unwrap_or_default(),
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
    reference_sequence: &fasta::record::Sequence,
    flags: Flags,
    alignment_start: Position,
    sequence: &Sequence,
    quality_scores: &QualityScores,
) -> io::Result<Vec<Feature>> {
    use noodles_sam::alignment::record::cigar::op::Kind;

    let mut features = Vec::new();

    let mut reference_position = alignment_start;
    let mut read_position = Position::MIN;

    for result in cigar.iter() {
        let op = result?;

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if op.len() == 1 {
                    let raw_reference_base = reference_sequence[reference_position];
                    let raw_read_base = sequence[read_position];

                    let quality_score = quality_scores[read_position];

                    if !flags.quality_scores_are_stored_as_array() {
                        features.push(Feature::QualityScore {
                            position: read_position,
                            quality_score,
                        });
                    }

                    if !raw_reference_base.eq_ignore_ascii_case(&raw_read_base) {
                        let feature = if let (Ok(reference_base), Ok(read_base)) = (
                            Base::try_from(raw_reference_base),
                            Base::try_from(raw_read_base),
                        ) {
                            Feature::Substitution {
                                position: read_position,
                                reference_base,
                                read_base,
                            }
                        } else {
                            Feature::ReadBase {
                                position: read_position,
                                base: raw_read_base,
                                quality_score,
                            }
                        };

                        features.push(feature);
                    }
                } else {
                    let read_end = read_position
                        .checked_add(op.len())
                        .expect("attempt to add with overflow");

                    if !flags.quality_scores_are_stored_as_array() {
                        let quality_scores = quality_scores[read_position..read_end].to_vec();

                        features.push(Feature::Scores {
                            position: read_position,
                            quality_scores,
                        });
                    }

                    let reference_end = reference_position
                        .checked_add(op.len())
                        .expect("attempt to add with overflow");

                    let reference_bases = &reference_sequence[reference_position..reference_end];
                    let read_bases = &sequence[read_position..read_end];

                    for (i, (&raw_reference_base, &raw_read_base)) in
                        reference_bases.iter().zip(read_bases).enumerate()
                    {
                        if !raw_reference_base.eq_ignore_ascii_case(&raw_read_base) {
                            let position = read_position
                                .checked_add(i)
                                .expect("attempt to add with overflow");

                            let feature = if let (Ok(reference_base), Ok(read_base)) = (
                                Base::try_from(raw_reference_base),
                                Base::try_from(raw_read_base),
                            ) {
                                Feature::Substitution {
                                    position,
                                    reference_base,
                                    read_base,
                                }
                            } else {
                                let quality_score = quality_scores[read_position];

                                Feature::ReadBase {
                                    position,
                                    base: raw_read_base,
                                    quality_score,
                                }
                            };

                            features.push(feature);
                        }
                    }
                }
            }
            Kind::Insertion => {
                if op.len() == 1 {
                    let base = sequence[read_position];
                    features.push(Feature::InsertBase {
                        position: read_position,
                        base,
                    });

                    if !flags.quality_scores_are_stored_as_array() {
                        let quality_score = quality_scores[read_position];

                        features.push(Feature::QualityScore {
                            position: read_position,
                            quality_score,
                        });
                    }
                } else {
                    let end = read_position
                        .checked_add(op.len())
                        .expect("attempt to add with overflow");

                    let bases = sequence[read_position..end].to_vec();
                    features.push(Feature::Insertion {
                        position: read_position,
                        bases,
                    });

                    if !flags.quality_scores_are_stored_as_array() {
                        let quality_scores = quality_scores[read_position..end].to_vec();

                        features.push(Feature::Scores {
                            position: read_position,
                            quality_scores,
                        });
                    }
                }
            }
            Kind::Deletion => features.push(Feature::Deletion {
                position: read_position,
                len: op.len(),
            }),
            Kind::Skip => features.push(Feature::ReferenceSkip {
                position: read_position,
                len: op.len(),
            }),
            Kind::SoftClip => {
                let end = read_position
                    .checked_add(op.len())
                    .expect("attempt to add with overflow");

                let bases = &sequence[read_position..end];

                features.push(Feature::SoftClip {
                    position: read_position,
                    bases: bases.to_vec(),
                });

                if !flags.quality_scores_are_stored_as_array() {
                    if bases.len() == 1 {
                        let quality_score = quality_scores[read_position];

                        features.push(Feature::QualityScore {
                            position: read_position,
                            quality_score,
                        });
                    } else {
                        let quality_scores = quality_scores[read_position..end].to_vec();

                        features.push(Feature::Scores {
                            position: read_position,
                            quality_scores,
                        });
                    }
                }
            }
            Kind::HardClip => features.push(Feature::HardClip {
                position: read_position,
                len: op.len(),
            }),
            Kind::Pad => features.push(Feature::Padding {
                position: read_position,
                len: op.len(),
            }),
        }

        if op.kind().consumes_reference() {
            reference_position = reference_position
                .checked_add(op.len())
                .expect("attempt to add with overflow");
        }

        if op.kind().consumes_read() {
            read_position = read_position
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
    use bstr::ByteSlice;
    use noodles_sam::alignment::{
        record::cigar::{Op, op::Kind},
        record_buf::Cigar,
    };

    use super::*;

    #[test]
    fn test_get_read_group_id() -> io::Result<()> {
        let read_group_name = b"rg0".as_bstr();

        let header = sam::Header::builder()
            .add_read_group(read_group_name, Default::default())
            .build();
        assert_eq!(get_read_group_id(&header, read_group_name)?, 0);

        let header = sam::Header::default();
        assert!(matches!(
            get_read_group_id(&header, read_group_name),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_cigar_to_features() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGT".to_vec());

        let flags = Flags::default();
        let alignment_start = Position::MIN;

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::QualityScore {
            position: Position::try_from(1)?,
            quality_score: 45,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"C");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::Scores {
            position: Position::try_from(1)?,
            quality_scores: vec![45, 35],
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AT");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::C,
                read_base: Base::T,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::InsertBase {
                position: Position::try_from(1)?,
                base: b'A',
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::QualityScore {
                position: Position::try_from(2)?,
                quality_score: 35,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::QualityScore {
                position: Position::try_from(3)?,
                quality_score: 43,
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Deletion {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::G,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Skip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::ReferenceSkip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::QualityScore {
                position: Position::try_from(2)?,
                quality_score: 35,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::QualityScore {
                position: Position::try_from(3)?,
                quality_score: 43,
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::HardClip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Pad, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Padding {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_cigar_to_features_with_quality_scores_stored_as_array()
    -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGT".to_vec());

        let flags = Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;
        let alignment_start = Position::MIN;

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"C");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::Substitution {
            position: Position::try_from(1)?,
            reference_base: Base::A,
            read_base: Base::C,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AT");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::Substitution {
            position: Position::try_from(2)?,
            reference_base: Base::C,
            read_base: Base::T,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::InsertBase {
                position: Position::try_from(1)?,
                base: b'A',
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Deletion {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::G,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Skip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::ReferenceSkip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::HardClip {
            position: Position::try_from(1)?,
            len: 1,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Pad, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            &reference_sequence,
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::Padding {
            position: Position::try_from(1)?,
            len: 1,
        }];
        assert_eq!(actual, expected);

        Ok(())
    }
}
