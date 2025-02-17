use std::cmp;

use crate::container::compression_header::preservation_map::SubstitutionMatrix;

/// A frequency table of reference-read base substitutions.
///
/// Each bin is the frequency of a substitution from the reference base to a read base. A same base
/// substitution is not a substitution; thus, the main diagonal will always be 0s.
///
/// ```text
///                         read base
///                     A   C   G   T   N
///                   +---+---+---+---+---+
///                 A | 0 |   |   |   |   |
///                   +---+---+---+---+---+
///                 C |   | 0 |   |   |   |
///                   +---+---+---+---+---+
/// reference base  G |   |   | 0 |   |   |
///                   +---+---+---+---+---+
///                 T |   |   |   | 0 |   |
///                   +---+---+---+---+---+
///                 N |   |   |   |   | 0 |
///                   +---+---+---+---+---+
/// ```
#[derive(Default)]
pub struct Frequencies([[u64; 5]; 5]);

impl Frequencies {
    pub fn hit(&mut self, reference_base: u8, read_base: u8) {
        let i = usize::from(encode(reference_base));
        let j = usize::from(encode(read_base));
        self.0[i][j] += 1;
    }

    fn row(&self, reference_base: u8) -> &[u64] {
        let i = usize::from(encode(reference_base));
        &self.0[i]
    }
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl From<u8> for Base {
    fn from(n: u8) -> Self {
        match n {
            b'A' => Self::A,
            b'C' => Self::C,
            b'G' => Self::G,
            b'T' => Self::T,
            _ => Self::N,
        }
    }
}

impl From<Frequencies> for SubstitutionMatrix {
    fn from(frequencies: Frequencies) -> Self {
        use crate::container::compression_header::preservation_map::substitution_matrix;

        const BASES: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

        let mut substitution_matrix = [[substitution_matrix::Base::N; 4]; 5];

        for reference_base in BASES {
            let mut base_frequencies: Vec<_> = BASES
                .into_iter()
                .zip(frequencies.row(reference_base).iter())
                .filter(|(read_base, _)| *read_base != reference_base)
                .collect();

            // § 10.6 "Mapped reads" (2021-10-15): "the substitutions for each reference base may
            // optionally be sorted by their frequencies, in descending order, with same-frequency
            // ties broken using the fixed order ACGTN."
            base_frequencies.sort_by_key(|(read_base, frequency)| {
                (cmp::Reverse(*frequency), Base::from(*read_base))
            });

            for (code, (read_base, _)) in base_frequencies.into_iter().enumerate() {
                let i = usize::from(encode(reference_base));
                // FIXME
                substitution_matrix[i][code] =
                    substitution_matrix::Base::try_from(read_base).unwrap();
            }
        }

        Self {
            substitutions: substitution_matrix,
        }
    }
}

fn encode(base: u8) -> u8 {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4, // N
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hit() {
        let mut frequencies = Frequencies::default();

        frequencies.hit(b'A', b'C');
        frequencies.hit(b'G', b'T');
        frequencies.hit(b'T', b'A');

        frequencies.hit(b'C', b'G');
        frequencies.hit(b'C', b'G');

        assert_eq!(
            frequencies.0,
            [
                [0, 1, 0, 0, 0],
                [0, 0, 2, 0, 0],
                [0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
            ]
        );
    }

    #[test]
    fn test_from_frequencies_for_substitution_matrix() {
        use crate::container::compression_header::preservation_map::substitution_matrix::Base;

        let frequencies = Frequencies([
            [0, 3, 8, 5, 0],
            [2, 0, 5, 5, 2],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        ]);

        let substitution_matrix = SubstitutionMatrix::from(frequencies);

        assert_eq!(
            substitution_matrix.substitutions,
            [
                [Base::G, Base::T, Base::C, Base::N],
                [Base::G, Base::T, Base::A, Base::N],
                [Base::A, Base::C, Base::T, Base::N],
                [Base::A, Base::C, Base::G, Base::N],
                [Base::A, Base::C, Base::G, Base::T],
            ]
        );
    }
}
