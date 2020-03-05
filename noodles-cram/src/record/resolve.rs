use std::convert::TryFrom;

use noodles_fasta as fasta;

use crate::{
    compression_header::preservation_map::substitution_matrix::Base,
    compression_header::SubstitutionMatrix, Feature,
};

pub fn resolve_bases(
    reference_sequence_record: &fasta::Record,
    substitution_matrix: &SubstitutionMatrix,
    features: &[Feature],
    alignment_start: i32,
    read_len: usize,
) -> Vec<u8> {
    let mut buf = vec![b'-'; read_len];

    let mut ref_pos = (alignment_start - 1) as usize;
    let mut read_pos = 0;

    let reference_sequence = reference_sequence_record.sequence();

    for feature in features {
        for _ in 0..(feature.position() - 1) {
            buf[read_pos] = reference_sequence[ref_pos];
            ref_pos += 1;
            read_pos += 1;
        }

        match feature {
            Feature::Substitution(_, code) => {
                let base = reference_sequence[ref_pos] as char;
                let reference_base =
                    Base::try_from(base).expect("invalid substitution matrix base");

                let read_base = substitution_matrix.get(reference_base, *code);
                buf[read_pos] = char::from(read_base) as u8;

                ref_pos += 1;
                read_pos += 1;
            }
            Feature::Deletion(_, len) => {
                ref_pos += *len as usize;

                for _ in 0..*len {
                    buf[read_pos] = reference_sequence[ref_pos];
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Feature::SoftClip(_, bases) => {
                let start = read_pos;
                let end = start + bases.len();
                buf.splice(start..end, bases.iter().cloned());

                read_pos += bases.len();
            }
            Feature::HardClip(..) => {}
            _ => todo!("resolve_bases: {:?}", feature),
        }
    }

    for base in buf.iter_mut().skip(read_pos) {
        *base = reference_sequence[ref_pos];
        ref_pos += 1;
    }

    buf
}
