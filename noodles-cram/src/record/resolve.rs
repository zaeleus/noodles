//! CRAM record field resolvers.

use std::io;

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, record::sequence::Base};

use crate::data_container::compression_header::SubstitutionMatrix;

use super::{
    feature::substitution::{self, Base as SubstitutionBase},
    Feature, Features,
};

pub(crate) fn resolve_bases(
    reference_sequence: Option<&fasta::record::Sequence>,
    substitution_matrix: &SubstitutionMatrix,
    features: &Features,
    alignment_start: Position,
    read_length: usize,
    buf: &mut sam::record::Sequence,
) -> io::Result<()> {
    buf.as_mut().fill(Base::N);
    buf.as_mut().resize(read_length, Base::N);

    let mut it = features.with_positions(alignment_start);

    let (mut last_reference_position, mut last_read_position) = it.positions();

    while let Some(((reference_position, read_position), feature)) = it.next() {
        if let Some(reference_sequence) = reference_sequence {
            let dst = &mut buf[last_read_position..read_position];
            let src = &reference_sequence[last_reference_position..reference_position];
            copy_from_raw_bases(dst, src)?;
        } else if read_position != last_read_position {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot resolve bases without reference sequence",
            ));
        }

        match feature {
            Feature::Bases(_, bases) => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Scores(..) => {}
            Feature::ReadBase(_, base, _) => buf[read_position] = *base,
            Feature::Substitution(_, substitution::Value::Code(code)) => {
                if let Some(reference_sequence) = reference_sequence {
                    let base = reference_sequence[reference_position];
                    let reference_base = SubstitutionBase::try_from(base).unwrap_or_default();
                    let read_base = substitution_matrix.get(reference_base, *code);
                    buf[read_position] = Base::from(read_base);
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "cannot resolve base substitution without reference sequence",
                    ));
                }
            }
            Feature::Substitution(_, substitution::Value::Bases(..)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "cannot resolve base substitution with bases",
                ))
            }
            Feature::Insertion(_, bases) => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Deletion(..) => {}
            Feature::InsertBase(_, base) => buf[read_position] = *base,
            Feature::QualityScore(..) => {}
            Feature::ReferenceSkip(..) => {}
            Feature::SoftClip(_, bases) => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Padding(..) => {}
            Feature::HardClip(..) => {}
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

        copy_from_raw_bases(dst, src)?;
    } else if buf[..last_read_position].len() != buf.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cannot resolve bases without reference sequence",
        ));
    }

    Ok(())
}

fn copy_from_bases(dst: &mut [Base], src: &[Base]) {
    for (&base, b) in src.iter().zip(dst.iter_mut()) {
        *b = base;
    }
}

fn copy_from_raw_bases(dst: &mut [Base], src: &[u8]) -> io::Result<()> {
    for (&raw_base, b) in src.iter().zip(dst.iter_mut()) {
        let base =
            Base::try_from(raw_base).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        *b = base;
    }

    Ok(())
}

/// Resolves the quality scores.
pub fn resolve_quality_scores(
    features: &[Feature],
    read_len: usize,
    quality_scores: &mut sam::record::QualityScores,
) {
    use sam::record::quality_scores::Score;

    quality_scores.as_mut().fill(Score::default());
    quality_scores.as_mut().resize(read_len, Score::default());

    for feature in features {
        let read_position = feature.position();

        match feature {
            Feature::Scores(_, scores) => {
                let end = read_position
                    .checked_add(scores.len())
                    .expect("attempt to add with overflow");

                quality_scores[read_position..end].copy_from_slice(scores);
            }
            Feature::ReadBase(_, _, score) => quality_scores[read_position] = *score,
            Feature::QualityScore(_, score) => quality_scores[read_position] = *score,
            _ => continue,
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;
    use sam::record::quality_scores::Score;

    use super::*;

    #[test]
    fn test_resolve_bases() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGTACGT".to_vec());
        let substitution_matrix = Default::default();
        let alignment_start = Position::try_from(1)?;

        let t = |features: &Features, expected: &sam::record::Sequence| {
            let mut actual = sam::record::Sequence::default();

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

        t(&Features::default(), &"ACGT".parse()?)?;
        t(
            &Features::from(vec![Feature::Bases(
                Position::try_from(1)?,
                vec![Base::T, Base::G],
            )]),
            &"TGGT".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::ReadBase(
                Position::try_from(2)?,
                Base::Y,
                Score::default(),
            )]),
            &"AYGT".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::Substitution(
                Position::try_from(2)?,
                substitution::Value::Code(1),
            )]),
            &"AGGT".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::Insertion(
                Position::try_from(2)?,
                vec![Base::G, Base::G],
            )]),
            &"AGGC".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::Deletion(Position::try_from(2)?, 2)]),
            &"ATAC".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::InsertBase(Position::try_from(2)?, Base::G)]),
            &"AGCG".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::ReferenceSkip(Position::try_from(2)?, 2)]),
            &"ATAC".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::SoftClip(
                Position::try_from(3)?,
                vec![Base::G, Base::G],
            )]),
            &"ACGG".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::Padding(Position::try_from(1)?, 2)]),
            &"ACGT".parse()?,
        )?;
        t(
            &Features::from(vec![Feature::HardClip(Position::try_from(1)?, 2)]),
            &"ACGT".parse()?,
        )?;

        let features = Features::from(vec![Feature::Substitution(
            Position::try_from(2)?,
            substitution::Value::Bases(SubstitutionBase::A, SubstitutionBase::C),
        )]);
        let mut actual = sam::record::Sequence::default();
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
        let features = Features::from(vec![Feature::Bases(
            Position::try_from(1)?,
            vec![Base::N, Base::N, Base::N, Base::N],
        )]);
        let alignment_start = Position::try_from(1)?;

        let mut actual = sam::record::Sequence::default();
        resolve_bases(
            None,
            &substitution_matrix,
            &features,
            alignment_start,
            4,
            &mut actual,
        )?;
        let expected = "NNNN".parse()?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_resolve_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        use sam::record::{quality_scores::Score, QualityScores};

        let features = [
            Feature::ReadBase(Position::try_from(1)?, Base::A, Score::try_from(5)?),
            Feature::QualityScore(Position::try_from(3)?, Score::try_from(8)?),
            Feature::Scores(
                Position::try_from(5)?,
                vec![Score::try_from(13)?, Score::try_from(21)?],
            ),
        ];

        let mut quality_scores = QualityScores::default();
        resolve_quality_scores(&features, 6, &mut quality_scores);
        let expected = QualityScores::try_from(vec![5, 0, 8, 0, 13, 21])?;
        assert_eq!(quality_scores, expected);

        Ok(())
    }
}
