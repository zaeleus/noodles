//! CRAM record field resolvers.

use std::io;

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, record::sequence::Base};

use crate::data_container::compression_header::SubstitutionMatrix;

use super::{Feature, Features};

pub(crate) fn resolve_bases(
    reference_sequence: Option<fasta::record::Sequence>,
    substitution_matrix: &SubstitutionMatrix,
    features: &Features,
    alignment_start: Position,
    read_length: usize,
) -> io::Result<sam::record::Sequence> {
    use crate::data_container::compression_header::preservation_map::substitution_matrix::Base as SubstitutionMatrixBase;

    let reference_sequence = reference_sequence.as_ref();

    let mut buf = sam::record::Sequence::from(vec![Base::N; read_length]);

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
            Feature::Substitution(_, code) => {
                if let Some(reference_sequence) = reference_sequence {
                    let base = reference_sequence[reference_position];
                    let reference_base = SubstitutionMatrixBase::try_from(base).unwrap_or_default();
                    let read_base = substitution_matrix.get(reference_base, *code);
                    buf[read_position] = Base::from(read_base);
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "cannot resolve base substitution without reference sequence",
                    ));
                }
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
    } else if usize::from(last_read_position) != buf.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cannot resolve bases without reference sequence",
        ));
    }

    Ok(buf)
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

/// Resolves the read features as CIGAR operations.
pub fn resolve_features(features: &Features, read_length: usize) -> sam::record::Cigar {
    use noodles_sam::record::cigar::{op::Kind, Op};

    let mut ops = Vec::new();
    let mut read_position = Position::MIN;

    for feature in features.iter() {
        if usize::from(feature.position()) > usize::from(read_position) {
            let len = usize::from(feature.position()) - usize::from(read_position);
            let op = Op::new(Kind::Match, len);
            ops.push(op);

            read_position = feature.position();
        }

        let (kind, len) = match feature {
            Feature::Substitution(..) => (Kind::Match, 1),
            Feature::Insertion(_, bases) => (Kind::Insertion, bases.len()),
            Feature::Deletion(_, len) => (Kind::Deletion, *len),
            Feature::InsertBase(..) => (Kind::Insertion, 1),
            Feature::ReferenceSkip(_, len) => (Kind::Skip, *len),
            Feature::SoftClip(_, bases) => (Kind::SoftClip, bases.len()),
            Feature::Padding(_, len) => (Kind::Pad, *len),
            Feature::HardClip(_, len) => (Kind::HardClip, *len),
            _ => continue,
        };

        let op = Op::new(kind, len);
        ops.push(op);

        if matches!(
            kind,
            Kind::Match
                | Kind::Insertion
                | Kind::SoftClip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
        ) {
            read_position = read_position
                .checked_add(len)
                .expect("attempt to add with overflow");
        }
    }

    if usize::from(read_position) <= read_length {
        let len = read_length - usize::from(read_position) + 1;
        let op = Op::new(Kind::Match, len);
        ops.push(op);
    }

    sam::record::Cigar::from(ops)
}

/// Resolves the quality scores.
pub fn resolve_quality_scores(features: &[Feature], read_len: usize) -> sam::record::QualityScores {
    use sam::record::quality_scores::Score;

    let mut scores = vec![Score::default(); read_len];

    for feature in features {
        let read_pos = usize::from(feature.position()) - 1;

        scores[read_pos] = match feature {
            Feature::ReadBase(_, _, quality_score) => *quality_score,
            Feature::QualityScore(_, quality_score) => *quality_score,
            _ => continue,
        }
    }

    sam::record::QualityScores::from(scores)
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
            let actual = resolve_bases(
                Some(reference_sequence.clone()),
                &substitution_matrix,
                features,
                alignment_start,
                4,
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
            &Features::from(vec![Feature::Substitution(Position::try_from(2)?, 1)]),
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

        Ok(())
    }

    #[test]
    fn test_resolve_features() -> Result<(), noodles_core::position::TryFromIntError> {
        use noodles_sam::record::{
            cigar::{op::Kind, Op},
            Cigar,
        };

        let features = Features::default();
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::Match, 4)])
        );

        let features = Features::from(vec![Feature::SoftClip(
            Position::try_from(1)?,
            vec![Base::A, Base::T],
        )]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 2)])
        );

        let features = Features::from(vec![Feature::SoftClip(
            Position::try_from(4)?,
            vec![Base::G],
        )]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::Match, 3), Op::new(Kind::SoftClip, 1)])
        );

        let features = Features::from(vec![Feature::HardClip(Position::try_from(1)?, 2)]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::HardClip, 2), Op::new(Kind::Match, 4)]),
        );

        // FIXME
        let features = Features::from(vec![
            Feature::SoftClip(Position::try_from(1)?, vec![Base::A]),
            Feature::Substitution(Position::try_from(3)?, 0),
        ]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![
                Op::new(Kind::SoftClip, 1),
                Op::new(Kind::Match, 1),
                Op::new(Kind::Match, 1),
                Op::new(Kind::Match, 1),
            ])
        );

        // FIXME
        let features = Features::from(vec![Feature::Substitution(Position::try_from(2)?, 0)]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Match, 1),
                Op::new(Kind::Match, 2)
            ])
        );

        Ok(())
    }

    #[test]
    fn test_resolve_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        use sam::record::{quality_scores::Score, QualityScores};

        let features = [
            Feature::ReadBase(Position::try_from(1)?, Base::A, Score::try_from(5)?),
            Feature::QualityScore(Position::try_from(3)?, Score::try_from(8)?),
        ];

        let actual = resolve_quality_scores(&features, 4);

        let expected = [5, 0, 8, 0]
            .into_iter()
            .map(Score::try_from)
            .collect::<Result<Vec<_>, _>>()
            .map(QualityScores::from)?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
