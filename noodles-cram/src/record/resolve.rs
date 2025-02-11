//! CRAM record field resolvers.

use std::io;

use noodles_core::Position;
use noodles_fasta as fasta;

use super::{
    feature::substitution::{self, Base as SubstitutionBase},
    Feature, Features, QualityScores, Sequence,
};
use crate::container::compression_header::SubstitutionMatrix;

pub(crate) fn resolve_bases(
    reference_sequence: Option<&fasta::record::Sequence>,
    substitution_matrix: &SubstitutionMatrix,
    features: &Features,
    alignment_start: Position,
    read_length: usize,
    buf: &mut Sequence,
) -> io::Result<()> {
    buf.as_mut().clear();
    buf.as_mut().resize(read_length, b'N');

    let mut it = features.with_positions(alignment_start);

    let (mut last_reference_position, mut last_read_position) = it.positions();

    while let Some(((reference_position, read_position), feature)) = it.next() {
        if let Some(reference_sequence) = reference_sequence {
            let dst = &mut buf[last_read_position..read_position];
            let src = &reference_sequence[last_reference_position..reference_position];
            copy_from_bases(dst, src);
        } else if read_position != last_read_position {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot resolve bases without reference sequence",
            ));
        }

        match feature {
            Feature::Bases { bases, .. } => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Scores { .. } => {}
            Feature::ReadBase { base, .. } => buf[read_position] = *base,
            Feature::Substitution { value, .. } => match value {
                substitution::Value::Code(code) => {
                    if let Some(reference_sequence) = reference_sequence {
                        let base = reference_sequence[reference_position];
                        let reference_base = SubstitutionBase::try_from(base).unwrap_or_default();
                        let read_base = substitution_matrix.get(reference_base, *code);
                        buf[read_position] = u8::from(read_base);
                    } else {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "cannot resolve base substitution without reference sequence",
                        ));
                    }
                }
                substitution::Value::Bases(..) => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "cannot resolve base substitution with bases",
                    ));
                }
            },
            Feature::Insertion { bases, .. } => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Deletion { .. } => {}
            Feature::InsertBase { base, .. } => buf[read_position] = *base,
            Feature::QualityScore { .. } => {}
            Feature::ReferenceSkip { .. } => {}
            Feature::SoftClip { bases, .. } => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Padding { .. } => {}
            Feature::HardClip { .. } => {}
        }

        let (next_reference_position, next_read_position) = it.positions();
        last_reference_position = next_reference_position;
        last_read_position = next_read_position;
    }

    if let Some(reference_sequence) = reference_sequence {
        let dst = &mut buf[last_read_position..];

        let end = last_reference_position
            .checked_add(dst.len())
            .expect("attempt to add with overflow");
        let src = &reference_sequence[last_reference_position..end];

        copy_from_bases(dst, src);
    } else if buf[..last_read_position].len() != buf.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cannot resolve bases without reference sequence",
        ));
    }

    Ok(())
}

fn copy_from_bases(dst: &mut [u8], src: &[u8]) {
    for (&base, b) in src.iter().zip(dst.iter_mut()) {
        *b = base;
    }
}

/// Resolves the quality scores.
pub fn resolve_quality_scores(
    features: &[Feature],
    read_len: usize,
    quality_scores: &mut QualityScores,
) {
    quality_scores.as_mut().clear();
    quality_scores.as_mut().resize(read_len, 0);

    for feature in features {
        let read_position = feature.position();

        match feature {
            Feature::Scores {
                quality_scores: scores,
                ..
            } => {
                let end = read_position
                    .checked_add(scores.len())
                    .expect("attempt to add with overflow");

                quality_scores[read_position..end].copy_from_slice(scores);
            }
            Feature::ReadBase {
                quality_score: score,
                ..
            } => quality_scores[read_position] = *score,
            Feature::QualityScore {
                quality_score: score,
                ..
            } => quality_scores[read_position] = *score,
            _ => continue,
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;

    #[test]
    fn test_resolve_bases() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGTACGT".to_vec());
        let substitution_matrix = Default::default();
        let alignment_start = Position::try_from(1)?;

        let t = |features: &Features, expected: &Sequence| {
            let mut actual = Sequence::default();

            resolve_bases(
                Some(&reference_sequence),
                &substitution_matrix,
                features,
                alignment_start,
                4,
                &mut actual,
            )?;

            assert_eq!(&actual, expected);

            Ok::<_, io::Error>(())
        };

        t(&Features::default(), &Sequence::from(b"ACGT"))?;
        t(
            &Features::from(vec![Feature::Bases {
                position: Position::try_from(1)?,
                bases: vec![b'T', b'G'],
            }]),
            &Sequence::from(b"TGGT"),
        )?;
        t(
            &Features::from(vec![Feature::ReadBase {
                position: Position::try_from(2)?,
                base: b'Y',
                quality_score: 0,
            }]),
            &Sequence::from(b"AYGT"),
        )?;
        t(
            &Features::from(vec![Feature::Substitution {
                position: Position::try_from(2)?,
                value: substitution::Value::Code(1),
            }]),
            &Sequence::from(b"AGGT"),
        )?;
        t(
            &Features::from(vec![Feature::Insertion {
                position: Position::try_from(2)?,
                bases: vec![b'G', b'G'],
            }]),
            &Sequence::from(b"AGGC"),
        )?;
        t(
            &Features::from(vec![Feature::Deletion {
                position: Position::try_from(2)?,
                len: 2,
            }]),
            &Sequence::from(b"ATAC"),
        )?;
        t(
            &Features::from(vec![Feature::InsertBase {
                position: Position::try_from(2)?,
                base: b'G',
            }]),
            &Sequence::from(b"AGCG"),
        )?;
        t(
            &Features::from(vec![Feature::ReferenceSkip {
                position: Position::try_from(2)?,
                len: 2,
            }]),
            &Sequence::from(b"ATAC"),
        )?;
        t(
            &Features::from(vec![Feature::SoftClip {
                position: Position::try_from(3)?,
                bases: vec![b'G', b'G'],
            }]),
            &Sequence::from(b"ACGG"),
        )?;
        t(
            &Features::from(vec![Feature::Padding {
                position: Position::try_from(1)?,
                len: 2,
            }]),
            &Sequence::from(b"ACGT"),
        )?;
        t(
            &Features::from(vec![Feature::HardClip {
                position: Position::try_from(1)?,
                len: 2,
            }]),
            &Sequence::from(b"ACGT"),
        )?;

        let features = Features::from(vec![Feature::Substitution {
            position: Position::try_from(2)?,
            value: substitution::Value::Bases(SubstitutionBase::A, SubstitutionBase::C),
        }]);
        let mut actual = Sequence::default();
        assert!(matches!(
            resolve_bases(
                Some(&reference_sequence),
                &substitution_matrix,
                &features,
                alignment_start,
                4,
                &mut actual,
            ),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_resolve_bases_without_a_reference_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let substitution_matrix = SubstitutionMatrix::default();
        let features = Features::from(vec![Feature::Bases {
            position: Position::try_from(1)?,
            bases: vec![b'N', b'N', b'N', b'N'],
        }]);
        let alignment_start = Position::try_from(1)?;

        let mut actual = Sequence::default();
        resolve_bases(
            None,
            &substitution_matrix,
            &features,
            alignment_start,
            4,
            &mut actual,
        )?;
        let expected = Sequence::from(b"NNNN");

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_resolve_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        let features = [
            Feature::ReadBase {
                position: Position::try_from(1)?,
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

        let mut quality_scores = QualityScores::default();
        resolve_quality_scores(&features, 6, &mut quality_scores);
        let expected = QualityScores::from(vec![5, 0, 8, 0, 13, 21]);
        assert_eq!(quality_scores, expected);

        Ok(())
    }
}
