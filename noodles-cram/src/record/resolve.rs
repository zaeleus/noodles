//! CRAM record field resolvers.

use std::io;

use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::data_container::compression_header::SubstitutionMatrix;

use super::{Feature, Features};

pub(crate) fn resolve_bases(
    reference_sequence: Option<fasta::record::Sequence>,
    substitution_matrix: &SubstitutionMatrix,
    features: &Features,
    alignment_start: sam::record::Position,
    read_length: usize,
) -> io::Result<Vec<u8>> {
    use crate::data_container::compression_header::preservation_map::substitution_matrix::Base;

    let reference_sequence = reference_sequence.as_ref().map(|rs| rs.as_ref());

    let mut buf = vec![b'-'; read_length];

    let mut it = features.with_positions(alignment_start);

    let (mut last_reference_position, mut last_read_position) = it.positions();
    last_reference_position -= 1;
    last_read_position -= 1;

    while let Some(((mut reference_position, mut read_position), feature)) = it.next() {
        reference_position -= 1;
        read_position -= 1;

        if let Some(reference_sequence) = reference_sequence {
            let dst = &mut buf[last_read_position..read_position];
            let src = &reference_sequence[last_reference_position..reference_position];
            dst.copy_from_slice(src);
        } else if read_position != last_read_position {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot resolve bases without reference sequence",
            ));
        }

        match feature {
            Feature::Bases(_, bases) => {
                let dst = &mut buf[read_position..read_position + bases.len()];
                dst.copy_from_slice(bases);
            }
            Feature::ReadBase(_, base, _) => {
                buf[read_position] = *base;
            }
            Feature::Substitution(_, code) => {
                if let Some(reference_sequence) = reference_sequence {
                    let base = char::from(reference_sequence[reference_position]);
                    let reference_base = Base::try_from(base).unwrap_or_default();
                    let read_base = substitution_matrix.get(reference_base, *code);
                    buf[read_position] = char::from(read_base) as u8;
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "cannot resolve base substitution without reference sequence",
                    ));
                }
            }
            Feature::Insertion(_, bases) => {
                let dst = &mut buf[read_position..read_position + bases.len()];
                dst.copy_from_slice(bases);
            }
            Feature::Deletion(..) => {}
            Feature::InsertBase(_, base) => {
                buf[read_position] = *base;
            }
            Feature::ReferenceSkip(..) => {}
            Feature::SoftClip(_, bases) => {
                let dst = &mut buf[read_position..read_position + bases.len()];
                dst.copy_from_slice(bases);
            }
            Feature::HardClip(..) => {}
            _ => todo!("resolve_bases: {:?}", feature),
        }

        let (next_reference_position, next_read_position) = it.positions();
        last_reference_position = next_reference_position - 1;
        last_read_position = next_read_position - 1;
    }

    if let Some(reference_sequence) = reference_sequence {
        let dst = &mut buf[last_read_position..];
        let src = &reference_sequence[last_reference_position..last_reference_position + dst.len()];
        dst.copy_from_slice(src);
    } else if last_read_position != buf.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "cannot resolve bases without reference sequence",
        ));
    }

    Ok(buf)
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
            Feature::Deletion(_, len) => (Kind::Deletion, *len),
            Feature::InsertBase(..) => (Kind::Insertion, 1),
            Feature::ReferenceSkip(_, len) => (Kind::Skip, *len),
            Feature::SoftClip(_, bases) => (Kind::SoftClip, bases.len() as i32),
            Feature::Padding(_, len) => (Kind::Pad, *len),
            Feature::HardClip(_, len) => (Kind::HardClip, *len),
            _ => continue,
        };

        let op = Op::new(kind, len as u32);
        ops.push(op);

        if matches!(
            kind,
            Kind::Match | Kind::Insertion | Kind::SoftClip | Kind::SeqMatch | Kind::SeqMismatch
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
pub fn resolve_quality_scores(features: &[Feature], read_len: usize) -> Vec<u8> {
    let mut quality_scores = vec![0; read_len];

    for feature in features {
        let read_pos = (feature.position() - 1) as usize;

        quality_scores[read_pos] = match feature {
            Feature::ReadBase(_, _, quality_score) => *quality_score,
            Feature::QualityScore(_, quality_score) => *quality_score,
            _ => continue,
        }
    }

    quality_scores
}

#[cfg(test)]
mod tests {
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
            &Features::from(vec![Feature::Bases(1, b"TG".to_vec())]),
            b"TGGT",
        )?;
        t(
            &Features::from(vec![Feature::ReadBase(2, b'Y', b'!')]),
            b"AYGT",
        )?;
        t(&Features::from(vec![Feature::Substitution(2, 1)]), b"AGGT")?;
        t(
            &Features::from(vec![Feature::Insertion(2, b"GG".to_vec())]),
            b"AGGC",
        )?;
        t(&Features::from(vec![Feature::Deletion(2, 2)]), b"ATAC")?;
        t(&Features::from(vec![Feature::InsertBase(2, b'G')]), b"AGCG")?;
        t(&Features::from(vec![Feature::ReferenceSkip(2, 2)]), b"ATAC")?;
        t(
            &Features::from(vec![Feature::SoftClip(3, b"GG".to_vec())]),
            b"ACGG",
        )?;
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

        let features = Features::from(vec![Feature::SoftClip(1, b"AT".to_vec())]);
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 2)])
        );

        let features = Features::from(vec![Feature::SoftClip(4, b"G".to_vec())]);
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
            Feature::SoftClip(1, b"A".to_vec()),
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
    fn test_resolve_quality_scores() {
        let features = [Feature::ReadBase(1, b'A', 5), Feature::QualityScore(3, 8)];
        assert_eq!(resolve_quality_scores(&features, 4), [5, 0, 8, 0]);
    }
}
