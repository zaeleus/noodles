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

    for feature in features {
        match feature {
            Feature::Substitution(position, code) => {
                let reference_sequence = reference_sequence_record.sequence();

                for _ in 0..(*position - 1) {
                    buf[read_pos] = reference_sequence[ref_pos];
                    ref_pos += 1;
                    read_pos += 1;
                }

                let base = reference_sequence[ref_pos] as char;
                let reference_base =
                    Base::try_from(base).expect("invalid substitution matrix base");

                let read_base = substitution_matrix.get(reference_base, *code);
                buf[read_pos] = char::from(read_base) as u8;

                ref_pos += 1;
                read_pos += 1;
            }
            Feature::SoftClip(position, bases) => {
                let reference_sequence = reference_sequence_record.sequence();

                for _ in 0..(*position - 1) {
                    buf[read_pos] = reference_sequence[ref_pos];
                    ref_pos += 1;
                    read_pos += 1;
                }

                let start = read_pos;
                let end = start + bases.len();
                buf.splice(start..end, bases.iter().cloned());

                read_pos += bases.len();
            }
            Feature::HardClip(position, _) => {
                let reference_sequence = reference_sequence_record.sequence();

                for _ in 0..(*position - 1) {
                    buf[read_pos] = reference_sequence[ref_pos];
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            _ => todo!("resolve_bases: {:?}", feature),
        }
    }

    buf
}
