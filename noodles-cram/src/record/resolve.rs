//! CRAM record field resolvers.

use std::io;

use noodles_fasta as fasta;
use noodles_sam::{self as sam, record::sequence::Base};

use crate::data_container::compression_header::SubstitutionMatrix;

use super::{Feature, Features};

pub(crate) fn resolve_bases(
    reference_sequence: Option<fasta::record::Sequence>,
    substitution_matrix: &SubstitutionMatrix,
    features: &Features,
    alignment_start: sam::record::Position,
    read_length: usize,
) -> io::Result<Vec<u8>> {
    use crate::data_container::compression_header::preservation_map::substitution_matrix::Base as SubstitutionMatrixBase;

    let reference_sequence = reference_sequence.as_ref();
    let mut buf = vec![b'-'; read_length];

    let mut it = features.with_positions(alignment_start);

    let (mut last_reference_position, mut last_read_position) = it.positions();
    last_read_position -= 1;

    while let Some(((reference_position, mut read_position), feature)) = it.next() {
        read_position -= 1;

        if let Some(reference_sequence) = reference_sequence {
            let dst = &mut buf[last_read_position..read_position];
            let src = reference_sequence
                .get(last_reference_position..reference_position)
                .unwrap();
            dst.copy_from_slice(src);
        } else if read_position != last_read_position {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot resolve bases without reference sequence",
            ));
        }

        match feature {
            Feature::Bases(_, bases) => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Scores(..) => {}
            Feature::ReadBase(_, base, _) => {
                buf[read_position] = u8::from(*base);
            }
            Feature::Substitution(_, code) => {
                if let Some(reference_sequence) = reference_sequence {
                    let base = reference_sequence.get(reference_position).copied().unwrap();
                    let reference_base = SubstitutionMatrixBase::try_from(base).unwrap_or_default();
                    let read_base = substitution_matrix.get(reference_base, *code);
                    buf[read_position] = u8::from(read_base);
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "cannot resolve base substitution without reference sequence",
                    ));
                }
            }
            Feature::Insertion(_, bases) => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Deletion(..) => {}
            Feature::InsertBase(_, base) => {
                buf[read_position] = u8::from(*base);
            }
            Feature::QualityScore(..) => {}
            Feature::ReferenceSkip(..) => {}
            Feature::SoftClip(_, bases) => copy_from_bases(&mut buf[read_position..], bases),
            Feature::Padding(..) => {}
            Feature::HardClip(..) => {}
        }

        let (next_reference_position, next_read_position) = it.positions();
        last_reference_position = next_reference_position;
        last_read_position = next_read_position - 1;
    }

    if let Some(reference_sequence) = reference_sequence {
        let dst = &mut buf[last_read_position..];

        let end = last_reference_position
            .checked_add(dst.len())
            .expect("attempt to add with overflow");
        let src = reference_sequence
            .get(last_reference_position..end)
            .unwrap();

        dst.copy_from_slice(src);
    } else if last_read_position != buf.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cannot resolve bases without reference sequence",
        ));
    }

    Ok(buf)
}

fn copy_from_bases(dst: &mut [u8], src: &[Base]) {
    for (&base, b) in src.iter().zip(dst.iter_mut()) {
        *b = u8::from(base);
    }
}

/// Resolves the read features as CIGAR operations.
pub fn resolve_features(features: &Features, read_len: i32) -> sam::record::Cigar {
    use noodles_sam::record::cigar::{op::Kind, Op};

    let mut ops = Vec::new();
    let mut i = 1;

    for feature in features.iter() {
        if feature.position() > i {
            let len = feature.position() - i;
            let op = Op::new(Kind::Match, len as u32);
            ops.push(op);

            i = feature.position();
        }

        let (kind, len) = match feature {
            Feature::Substitution(..) => (Kind::Match, 1),
            Feature::Insertion(_, bases) => (Kind::Insertion, bases.len() as i32),
            Feature::Deletion(_, len) => (Kind::Deletion, *len as i32),
            Feature::InsertBase(..) => (Kind::Insertion, 1),
            Feature::ReferenceSkip(_, len) => (Kind::Skip, *len as i32),
            Feature::SoftClip(_, bases) => (Kind::SoftClip, bases.len() as i32),
            Feature::Padding(_, len) => (Kind::Pad, *len as i32),
            Feature::HardClip(_, len) => (Kind::HardClip, *len as i32),
            _ => continue,
        };

        let op = Op::new(kind, len as u32);
        ops.push(op);

        if matches!(
            kind,
            Kind::Match
                | Kind::Insertion
                | Kind::SoftClip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
        ) {
            i += len;
        }
    }

    if i <= read_len {
        let len = read_len - i + 1;
        let op = Op::new(Kind::Match, len as u32);
        ops.push(op);
    }

    sam::record::Cigar::from(ops)
}

/// Resolves the quality scores.
pub fn resolve_quality_scores(features: &[Feature], read_len: usize) -> sam::record::QualityScores {
    use sam::record::quality_scores::Score;

    let mut scores = vec![Score::default(); read_len];

    for feature in features {
        let read_pos = (feature.position() - 1) as usize;

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
    use sam::record::quality_scores::Score;

    use super::*;

    #[test]
    fn test_resolve_bases() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGTACGT".to_vec());
        let substitution_matrix = Default::default();
        let alignment_start = sam::record::Position::try_from(1)?;

        let t = |features: &Features, expected: &[u8]| {
            let actual = resolve_bases(
                Some(reference_sequence.clone()),
                &substitution_matrix,
                features,
                alignment_start,
                4,
            )?;

            assert_eq!(actual, expected);

            Ok::<_, io::Error>(())
        };

        t(&Features::default(), b"ACGT")?;
        t(
            &Features::from(vec![Feature::Bases(1, vec![Base::T, Base::G])]),
            b"TGGT",
        )?;
        t(
            &Features::from(vec![Feature::ReadBase(2, Base::Y, Score::default())]),
            b"AYGT",
        )?;
        t(&Features::from(vec![Feature::Substitution(2, 1)]), b"AGGT")?;
        t(
            &Features::from(vec![Feature::Insertion(2, vec![Base::G, Base::G])]),
            b"AGGC",
        )?;
        t(&Features::from(vec![Feature::Deletion(2, 2)]), b"ATAC")?;
        t(
            &Features::from(vec![Feature::InsertBase(2, Base::G)]),
            b"AGCG",
        )?;
        t(&Features::from(vec![Feature::ReferenceSkip(2, 2)]), b"ATAC")?;
        t(
            &Features::from(vec![Feature::SoftClip(3, vec![Base::G, Base::G])]),
            b"ACGG",
        )?;
        t(&Features::from(vec![Feature::Padding(1, 2)]), b"ACGT")?;
        t(&Features::from(vec![Feature::HardClip(1, 2)]), b"ACGT")?;

        Ok(())
    }

    #[test]
    fn test_resolve_features() {
        use noodles_sam::record::{
            cigar::{op::Kind, Op},
            Cigar,
        };

        let features = Features::default();
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::Match, 4)])
        );

        let features = Features::from(vec![Feature::SoftClip(1, vec![Base::A, Base::T])]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 2)])
        );

        let features = Features::from(vec![Feature::SoftClip(4, vec![Base::G])]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::Match, 3), Op::new(Kind::SoftClip, 1)])
        );

        let features = Features::from(vec![Feature::HardClip(1, 2)]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::HardClip, 2), Op::new(Kind::Match, 4)]),
        );

        // FIXME
        let features = Features::from(vec![
            Feature::SoftClip(1, vec![Base::A]),
            Feature::Substitution(3, 0),
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
        let features = Features::from(vec![Feature::Substitution(2, 0)]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Match, 1),
                Op::new(Kind::Match, 2)
            ])
        );
    }

    #[test]
    fn test_resolve_quality_scores(
    ) -> Result<(), sam::record::quality_scores::score::TryFromUByteError> {
        use sam::record::{quality_scores::Score, QualityScores};

        let features = [
            Feature::ReadBase(1, Base::A, Score::try_from(5)?),
            Feature::QualityScore(3, Score::try_from(8)?),
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
