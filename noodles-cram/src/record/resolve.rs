//! CRAM record field resolvers.

use noodles_fasta as fasta;
use noodles_sam::record::Cigar;

use super::Feature;
use crate::data_container::CompressionHeader;

/// Resolves the read bases.
pub fn resolve_bases(
    reference_sequence_record: &fasta::Record,
    compression_header: &CompressionHeader,
    features: &[Feature],
    alignment_start: i32,
    read_len: usize,
) -> Vec<u8> {
    let reference_sequence = reference_sequence_record.sequence();
    let substitution_matrix = compression_header.preservation_map().substitution_matrix();

    internal::resolve_bases(
        reference_sequence,
        substitution_matrix,
        features,
        alignment_start,
        read_len,
    )
}

mod internal {
    use crate::data_container::compression_header::preservation_map::SubstitutionMatrix;

    use super::*;

    pub(super) fn resolve_bases(
        reference_sequence: &fasta::record::Sequence,
        substitution_matrix: &SubstitutionMatrix,
        features: &[Feature],
        alignment_start: i32,
        read_length: usize,
    ) -> Vec<u8> {
        use crate::data_container::compression_header::preservation_map::substitution_matrix::Base;

        let raw_reference_sequence = reference_sequence.as_ref();

        let mut buf = vec![b'-'; read_length];

        let mut ref_pos = (alignment_start - 1) as usize;
        let mut read_pos = 0;

        for feature in features {
            let feature_pos = feature.position() as usize;

            while read_pos < feature_pos - 1 {
                buf[read_pos] = raw_reference_sequence[ref_pos];
                ref_pos += 1;
                read_pos += 1;
            }

            match feature {
                Feature::Substitution(_, code) => {
                    let base = raw_reference_sequence[ref_pos] as char;
                    let reference_base = Base::try_from(base).unwrap_or_default();

                    let read_base = substitution_matrix.get(reference_base, *code);
                    buf[read_pos] = char::from(read_base) as u8;

                    ref_pos += 1;
                    read_pos += 1;
                }
                Feature::Insertion(_, bases) => {
                    for &base in bases {
                        buf[read_pos] = base;
                        read_pos += 1;
                    }
                }
                Feature::Deletion(_, len) => {
                    ref_pos += *len as usize;
                }
                Feature::InsertBase(_, base) => {
                    buf[read_pos] = *base;
                    read_pos += 1;
                }
                Feature::ReferenceSkip(_, len) => {
                    ref_pos += *len as usize;
                }
                Feature::SoftClip(_, bases) => {
                    for &base in bases {
                        buf[read_pos] = base;
                        read_pos += 1;
                    }
                }
                Feature::HardClip(..) => {}
                _ => todo!("resolve_bases: {:?}", feature),
            }
        }

        for base in buf.iter_mut().skip(read_pos) {
            *base = raw_reference_sequence[ref_pos];
            ref_pos += 1;
        }

        buf
    }
}

/// Resolves the read features as CIGAR operations.
pub fn resolve_features(features: &[Feature], read_len: i32) -> Cigar {
    use noodles_sam::record::cigar::{op::Kind, Op};

    let mut ops = Vec::new();
    let mut i = 1;

    for feature in features {
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

    Cigar::from(ops)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolve_bases() {
        let reference_sequence = fasta::record::Sequence::from(b"ACGTACGT".to_vec());
        let substitution_matrix = Default::default();

        let t = |features: &[Feature], expected: &[u8]| {
            let actual =
                internal::resolve_bases(&reference_sequence, &substitution_matrix, features, 1, 4);
            assert_eq!(actual, expected);
        };

        t(&[], b"ACGT");
        t(&[Feature::Substitution(2, 1)], b"AGGT");
        t(&[Feature::Insertion(2, b"GG".to_vec())], b"AGGC");
        t(&[Feature::Deletion(2, 2)], b"ATAC");
        t(&[Feature::InsertBase(2, b'G')], b"AGCG");
        t(&[Feature::ReferenceSkip(2, 2)], b"ATAC");
        t(&[Feature::SoftClip(3, b"GG".to_vec())], b"ACGG");
        t(&[Feature::HardClip(1, 2)], b"ACGT");
    }

    #[test]
    fn test_resolve_features() {
        use noodles_sam::record::cigar::{op::Kind, Op};

        let features = [];
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::Match, 4)])
        );

        let features = [Feature::SoftClip(1, b"AT".to_vec())];
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 2)])
        );

        let features = [Feature::SoftClip(4, b"G".to_vec())];
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::Match, 3), Op::new(Kind::SoftClip, 1)])
        );

        let features = [Feature::HardClip(1, 2)];
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![Op::new(Kind::HardClip, 2), Op::new(Kind::Match, 4)]),
        );

        // FIXME
        let features = [
            Feature::SoftClip(1, b"A".to_vec()),
            Feature::Substitution(3, 0),
        ];
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
        let features = [Feature::Substitution(2, 0)];
        assert_eq!(
            resolve_features(&features, 4),
            Cigar::from(vec![
                Op::new(Kind::Match, 1),
                Op::new(Kind::Match, 1),
                Op::new(Kind::Match, 2)
            ])
        );
    }
}
