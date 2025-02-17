mod frequencies;

use std::io::{self, Write};

use self::frequencies::Frequencies;
use crate::{
    container::compression_header::preservation_map::{
        substitution_matrix::Base, SubstitutionMatrix,
    },
    io::writer::{record::Feature, Record},
};

pub(super) fn write_substitution_matrix<W>(
    writer: &mut W,
    substitution_matrix: &SubstitutionMatrix,
) -> io::Result<()>
where
    W: Write,
{
    let buf = encode(substitution_matrix);
    writer.write_all(&buf)
}

fn encode(substitution_matrix: &SubstitutionMatrix) -> [u8; 5] {
    let mut buf = [0; 5];
    let mut index_bases = [(0, Base::N); 4];

    for (bases, codes) in substitution_matrix.substitutions.iter().zip(&mut buf) {
        for ((i, base), index_base) in bases.iter().enumerate().zip(&mut index_bases) {
            *index_base = (i, *base);
        }

        index_bases.sort_by_key(|(_, base)| *base);

        for (i, _) in &index_bases {
            *codes <<= 2;
            // SAFETY: `i < 4`.
            *codes |= *i as u8;
        }
    }

    buf
}

pub(super) fn build_substitution_matrix(records: &[Record]) -> SubstitutionMatrix {
    let mut frequencies = Frequencies::default();

    for record in records {
        for feature in &record.features {
            if let Feature::Substitution {
                reference_base,
                read_base,
                ..
            } = feature
            {
                frequencies.hit(*reference_base, *read_base);
            }
        }
    }

    frequencies.into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() {
        // § 10.6.4 "Mapped reads: Substitution Matrix Format" (2024-09-04)
        let substitution_matrix = SubstitutionMatrix {
            substitutions: [
                [Base::T, Base::C, Base::G, Base::N], // A
                [Base::G, Base::A, Base::T, Base::N], // C
                [Base::C, Base::T, Base::A, Base::N], // G
                [Base::A, Base::G, Base::C, Base::N], // T
                [Base::A, Base::C, Base::G, Base::T], // N
            ],
        };

        assert_eq!(encode(&substitution_matrix), [0x63, 0x4b, 0x87, 0x27, 0x1b]);
    }
}
