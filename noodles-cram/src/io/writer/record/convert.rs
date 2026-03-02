use std::io;

use bstr::BStr;
use noodles_core::Position;
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
    record::{self, Flags, MateFlags},
};

impl Record {
    pub fn try_from_alignment_record_with_options(
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
        strip_md_nm: bool,
    ) -> io::Result<Self> {
        let bam_flags = record.flags()?;
        let mut cram_flags = Flags::default();

        let sequence = Sequence::from(record.sequence().iter().collect::<Vec<_>>());

        let quality_scores = if record.quality_scores().is_empty() {
            cram_flags.insert(Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY);
            vec![0xff; record.sequence().len()].into()
        } else {
            if bam_flags.is_unmapped() {
                cram_flags.insert(Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY);
            }

            record.quality_scores().iter().collect::<io::Result<_>>()?
        };

        let features = cigar_to_features(
            record.cigar().as_ref(),
            cram_flags,
            &sequence,
            &quality_scores,
        )?;

        let data = record.data();
        let (data_buf, read_group_name) = get_filtered_data(data.as_ref(), strip_md_nm)?;

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
            template_length: i64::from(record.template_length()?),
            mate_distance: None,
            data: data_buf,
            features,
            mapping_quality: record.mapping_quality().transpose()?,
            sequence: sequence.into(),
            quality_scores: quality_scores.into(),
        })
    }

    /// Creates a writer Record directly from a CRAM reader Record, avoiding the
    /// CIGAR round-trip that occurs when going through the `sam::alignment::Record` trait.
    pub fn from_cram_record(
        cram_record: &crate::Record<'_>,
        strip_md_nm: bool,
    ) -> io::Result<Self> {
        let features = cram_record
            .features
            .iter()
            .map(|f| convert_feature(f, &cram_record.substitution_matrix, cram_record))
            .collect::<io::Result<Vec<_>>>()?;

        let data_buf = convert_cram_data(&cram_record.data, strip_md_nm)?;

        Ok(Self {
            bam_flags: cram_record.bam_flags,
            cram_flags: cram_record.cram_flags,
            reference_sequence_id: cram_record.reference_sequence_id,
            read_length: cram_record.read_length,
            alignment_start: cram_record.alignment_start,
            read_group_id: cram_record.read_group_id,
            name: cram_record.name.as_deref().map(|s| s.into()),
            mate_flags: MateFlags::default(),
            mate_reference_sequence_id: cram_record.mate_reference_sequence_id,
            mate_alignment_start: cram_record.mate_alignment_start,
            template_length: cram_record.template_length,
            mate_distance: None,
            data: data_buf,
            features,
            mapping_quality: cram_record.mapping_quality,
            sequence: cram_record.sequence.to_vec(),
            quality_scores: cram_record.quality_scores.to_vec(),
        })
    }
}

fn convert_feature(
    feature: &record::Feature<'_>,
    substitution_matrix: &crate::container::compression_header::preservation_map::SubstitutionMatrix,
    record: &crate::Record<'_>,
) -> io::Result<Feature> {
    Ok(match feature {
        record::Feature::Bases { position, bases } => Feature::Bases {
            position: *position,
            bases: bases.to_vec(),
        },
        record::Feature::Scores {
            position,
            quality_scores,
        } => Feature::Scores {
            position: *position,
            quality_scores: quality_scores.to_vec(),
        },
        record::Feature::ReadBase {
            position,
            base,
            quality_score,
        } => Feature::ReadBase {
            position: *position,
            base: *base,
            quality_score: *quality_score,
        },
        record::Feature::Substitution { position, code } => {
            // Resolve substitution code to actual bases using the substitution matrix
            // and reference sequence
            let ref_base_u8 = get_reference_base_at(record, *position)?;
            let ref_base = Base::try_from(ref_base_u8).unwrap_or(Base::N);
            let read_base = substitution_matrix.get(ref_base, *code);

            Feature::Substitution {
                position: *position,
                reference_base: ref_base,
                read_base,
            }
        }
        record::Feature::Insertion { position, bases } => Feature::Insertion {
            position: *position,
            bases: bases.to_vec(),
        },
        record::Feature::Deletion { position, len } => Feature::Deletion {
            position: *position,
            len: *len,
        },
        record::Feature::InsertBase { position, base } => Feature::InsertBase {
            position: *position,
            base: *base,
        },
        record::Feature::QualityScore {
            position,
            quality_score,
        } => Feature::QualityScore {
            position: *position,
            quality_score: *quality_score,
        },
        record::Feature::ReferenceSkip { position, len } => Feature::ReferenceSkip {
            position: *position,
            len: *len,
        },
        record::Feature::SoftClip { position, bases } => Feature::SoftClip {
            position: *position,
            bases: bases.to_vec(),
        },
        record::Feature::Padding { position, len } => Feature::Padding {
            position: *position,
            len: *len,
        },
        record::Feature::HardClip { position, len } => Feature::HardClip {
            position: *position,
            len: *len,
        },
    })
}

fn get_reference_base_at(record: &crate::Record<'_>, read_position: Position) -> io::Result<u8> {
    use crate::io::reader::container::slice::ReferenceSequence;

    let alignment_start = match record.alignment_start {
        Some(pos) => usize::from(pos),
        None => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "reference base requested but alignment start is missing",
            ));
        }
    };

    // Calculate the reference offset by walking features up to this position
    let mut ref_offset: usize = 0;
    let mut read_pos: usize = 1;
    let feat_pos = usize::from(read_position);

    for f in &record.features {
        let fp = usize::from(f.position());
        if fp >= feat_pos {
            break;
        }
        match f {
            record::Feature::Insertion { bases, .. } => {
                let matches_before = fp - read_pos;
                ref_offset += matches_before;
                read_pos = fp + bases.len();
                // Insertion doesn't consume reference
            }
            record::Feature::InsertBase { .. } => {
                let matches_before = fp - read_pos;
                ref_offset += matches_before;
                read_pos = fp + 1;
            }
            record::Feature::Deletion { len, .. } => {
                let matches_before = fp - read_pos;
                ref_offset += matches_before;
                ref_offset += len;
                read_pos = fp;
            }
            record::Feature::ReferenceSkip { len, .. } => {
                let matches_before = fp - read_pos;
                ref_offset += matches_before;
                ref_offset += len;
                read_pos = fp;
            }
            record::Feature::SoftClip { bases, .. } => {
                let matches_before = fp - read_pos;
                ref_offset += matches_before;
                read_pos = fp + bases.len();
            }
            record::Feature::HardClip { .. }
            | record::Feature::Padding { .. }
            | record::Feature::Scores { .. }
            | record::Feature::QualityScore { .. } => {
                // These don't consume read or reference positions
            }
            _ => {
                // For features that consume both read and ref equally
                let matches_before = fp - read_pos;
                ref_offset += matches_before;
                match f {
                    record::Feature::Bases { bases, .. } => {
                        ref_offset += bases.len();
                        read_pos = fp + bases.len();
                    }
                    record::Feature::ReadBase { .. } | record::Feature::Substitution { .. } => {
                        ref_offset += 1;
                        read_pos = fp + 1;
                    }
                    _ => unreachable!("all feature variants should be handled explicitly"),
                }
            }
        }
    }

    // Add remaining matches between last feature and current position
    ref_offset += feat_pos - read_pos;

    // ref_index is 0-based genomic position
    let ref_index = alignment_start - 1 + ref_offset;

    match &record.reference_sequence {
        Some(ReferenceSequence::External { sequence }) => {
            let seq: &[u8] = (**sequence).as_ref();
            seq.get(ref_index).copied().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "reference base index {ref_index} out of bounds (len: {})",
                        seq.len()
                    ),
                )
            })
        }
        Some(ReferenceSequence::Embedded {
            reference_start,
            sequence,
        }) => {
            let start_0based = usize::from(*reference_start) - 1;
            let adjusted = ref_index.checked_sub(start_0based).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "reference index {ref_index} is before embedded sequence start {start_0based}"
                    ),
                )
            })?;
            sequence.get(adjusted).copied().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "embedded reference base index {adjusted} out of bounds (len: {})",
                        sequence.len()
                    ),
                )
            })
        }
        None => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "reference base requested for Substitution feature but no reference sequence is available",
        )),
    }
}

fn convert_cram_data(
    data: &[(Tag, ValueBuf)],
    strip_md_nm: bool,
) -> io::Result<Vec<(Tag, ValueBuf)>> {
    if !strip_md_nm {
        return Ok(data.to_vec());
    }

    let mut data_buf = Vec::with_capacity(data.len());

    for (tag, value) in data {
        if *tag == Tag::MISMATCHED_POSITIONS || *tag == Tag::EDIT_DISTANCE {
            continue;
        }

        data_buf.push((*tag, value.clone()));
    }

    Ok(data_buf)
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
    strip_md_nm: bool,
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

        if strip_md_nm && (tag == Tag::MISMATCHED_POSITIONS || tag == Tag::EDIT_DISTANCE) {
            continue;
        }

        data_buf.push((tag, value.try_into()?));
    }

    Ok((data_buf, read_group_name))
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::{
        record::cigar::{Op, op::Kind},
        record_buf::Cigar,
    };

    use super::*;

    #[test]
    fn test_cigar_to_features() -> Result<(), Box<dyn std::error::Error>> {
        let flags = Flags::default();

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![Feature::ReadBase {
            position: Position::try_from(1)?,
            base: b'A',
            quality_score: 45,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
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
        let quality_scores = [45, 35].into_iter().collect();
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
        let quality_scores = [45, 35, 43].into_iter().collect();
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
        let quality_scores = [45, 35].into_iter().collect();
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
        let quality_scores = [45].into_iter().collect();
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
        let quality_scores = [45, 35].into_iter().collect();
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
        let quality_scores = [45, 35, 43].into_iter().collect();
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
        let quality_scores = [45].into_iter().collect();
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
        let quality_scores = [45].into_iter().collect();
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
    fn test_cigar_to_features_with_quality_scores_stored_as_array()
    -> Result<(), Box<dyn std::error::Error>> {
        let flags = Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(&cigar, flags, &sequence, &quality_scores)?;
        let expected = vec![Feature::ReadBase {
            position: Position::try_from(1)?,
            base: b'A',
            quality_score: 45,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
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
        let quality_scores = [45, 35].into_iter().collect();
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
        let quality_scores = [45, 35, 43].into_iter().collect();
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
        let quality_scores = [45, 35].into_iter().collect();
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
        let quality_scores = [45].into_iter().collect();
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
        let quality_scores = [45, 35].into_iter().collect();
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
        let quality_scores = [45, 35, 43].into_iter().collect();
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
        let quality_scores = [45].into_iter().collect();
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
        let quality_scores = [45].into_iter().collect();
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
    fn test_from_cram_record() -> Result<(), Box<dyn std::error::Error>> {
        use std::borrow::Cow;

        use noodles_fasta::record::Sequence as FastaSequence;
        use noodles_sam::alignment::record::MappingQuality;

        use crate::{
            container::compression_header::preservation_map::SubstitutionMatrix,
            io::reader::container::slice::ReferenceSequence,
            record::{Feature as CramFeature, Flags as CramFlags},
        };

        // Build a minimal reference sequence: "ACGTACGT"
        let ref_seq = FastaSequence::from(b"ACGTACGT".to_vec());

        let mut cram_record = crate::Record::default();
        cram_record.bam_flags = sam::alignment::record::Flags::empty();
        cram_record.cram_flags = CramFlags::default();
        cram_record.reference_sequence_id = Some(0);
        cram_record.read_length = 4;
        cram_record.alignment_start = Position::new(1);
        cram_record.name = Some(Cow::Borrowed(b"read1"));
        cram_record.sequence = Cow::Borrowed(b"ACGT");
        cram_record.quality_scores = Cow::Borrowed(&[40, 35, 30, 25]);
        cram_record.mapping_quality = MappingQuality::new(30);
        cram_record.reference_sequence = Some(ReferenceSequence::External {
            sequence: ref_seq.into(),
        });
        cram_record.substitution_matrix = SubstitutionMatrix::default();

        // Add a deletion feature at position 3
        cram_record.features = vec![CramFeature::Deletion {
            position: Position::try_from(3)?,
            len: 2,
        }];

        let writer_record = super::Record::from_cram_record(&cram_record, false)?;

        assert_eq!(writer_record.reference_sequence_id, Some(0));
        assert_eq!(writer_record.read_length, 4);
        assert_eq!(writer_record.alignment_start, Position::new(1));
        assert!(writer_record.name.as_ref().map(|n| n.as_ref()) == Some(b"read1".as_slice()));
        assert_eq!(writer_record.sequence, b"ACGT");
        assert_eq!(writer_record.quality_scores, vec![40, 35, 30, 25]);
        assert_eq!(writer_record.mapping_quality, MappingQuality::new(30));
        assert_eq!(writer_record.features.len(), 1);
        assert_eq!(
            writer_record.features[0],
            Feature::Deletion {
                position: Position::try_from(3)?,
                len: 2,
            }
        );

        Ok(())
    }

    #[test]
    fn test_from_cram_record_with_substitution() -> Result<(), Box<dyn std::error::Error>> {
        use std::borrow::Cow;

        use noodles_fasta::record::Sequence as FastaSequence;

        use crate::{
            container::compression_header::preservation_map::{
                SubstitutionMatrix, substitution_matrix::Base,
            },
            io::reader::container::slice::ReferenceSequence,
            record::{Feature as CramFeature, Flags as CramFlags},
        };

        // Reference: "ACGTACGT", read has substitution at position 2 (C -> T)
        let ref_seq = FastaSequence::from(b"ACGTACGT".to_vec());
        let sub_matrix = SubstitutionMatrix::default();
        // With default matrix, for ref base C, code 0b10 => T
        let sub_code = sub_matrix.find(Base::C, Base::T);

        let mut cram_record = crate::Record::default();
        cram_record.bam_flags = sam::alignment::record::Flags::empty();
        cram_record.cram_flags = CramFlags::default();
        cram_record.reference_sequence_id = Some(0);
        cram_record.read_length = 4;
        cram_record.alignment_start = Position::new(1);
        cram_record.sequence = Cow::Borrowed(b"ATGT");
        cram_record.quality_scores = Cow::Borrowed(&[40, 35, 30, 25]);
        cram_record.reference_sequence = Some(ReferenceSequence::External {
            sequence: ref_seq.into(),
        });
        cram_record.substitution_matrix = sub_matrix;

        cram_record.features = vec![CramFeature::Substitution {
            position: Position::try_from(2)?,
            code: sub_code,
        }];

        let writer_record = super::Record::from_cram_record(&cram_record, false)?;

        assert_eq!(writer_record.features.len(), 1);
        assert_eq!(
            writer_record.features[0],
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::C,
                read_base: Base::T,
            }
        );

        Ok(())
    }

    #[test]
    fn test_from_cram_record_strip_md_nm() -> Result<(), Box<dyn std::error::Error>> {
        use std::borrow::Cow;

        use noodles_sam::alignment::record_buf::data::field::Value as ValueBuf;

        use crate::record::Flags as CramFlags;

        let mut cram_record = crate::Record::default();
        cram_record.bam_flags = sam::alignment::record::Flags::UNMAPPED;
        cram_record.cram_flags = CramFlags::default();
        cram_record.sequence = Cow::Borrowed(b"ACGT");
        cram_record.quality_scores = Cow::Borrowed(&[40, 35, 30, 25]);

        // Add MD and NM tags along with a normal tag
        cram_record.data = vec![
            (Tag::MISMATCHED_POSITIONS, ValueBuf::from("4")),
            (Tag::EDIT_DISTANCE, ValueBuf::from(0u32)),
            (Tag::ALIGNMENT_SCORE, ValueBuf::from(100i32)),
        ];

        // Without stripping
        let writer_record = super::Record::from_cram_record(&cram_record, false)?;
        assert_eq!(writer_record.data.len(), 3);

        // With stripping
        let writer_record = super::Record::from_cram_record(&cram_record, true)?;
        assert_eq!(writer_record.data.len(), 1);
        assert_eq!(writer_record.data[0].0, Tag::ALIGNMENT_SCORE);

        Ok(())
    }
}
