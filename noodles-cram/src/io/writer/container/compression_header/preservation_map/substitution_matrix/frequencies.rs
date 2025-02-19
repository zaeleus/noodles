use std::cmp;

use crate::container::compression_header::preservation_map::{
    substitution_matrix::Base, SubstitutionMatrix,
};

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
    pub fn hit(&mut self, reference_base: Base, read_base: Base) {
        let i = reference_base as usize;
        let j = read_base as usize;
        self.0[i][j] += 1;
    }

    fn row(&self, reference_base: Base) -> &[u64; 5] {
        let i = reference_base as usize;
        &self.0[i]
    }
}

impl From<Frequencies> for SubstitutionMatrix {
    fn from(frequencies: Frequencies) -> Self {
        const BASES: [Base; 5] = [Base::A, Base::C, Base::G, Base::T, Base::N];

        let mut substitution_matrix = [[Base::N; 4]; 5];

        for reference_base in BASES {
            let mut base_frequencies: Vec<_> = BASES
                .into_iter()
                .zip(frequencies.row(reference_base).iter())
                .filter(|(read_base, _)| *read_base != reference_base)
                .collect();

            // § 10.6 "Mapped reads" (2021-10-15): "the substitutions for each reference base may
            // optionally be sorted by their frequencies, in descending order, with same-frequency
            // ties broken using the fixed order ACGTN."
            base_frequencies
                .sort_by_key(|(read_base, frequency)| (cmp::Reverse(*frequency), *read_base));

            for (code, (read_base, _)) in base_frequencies.into_iter().enumerate() {
                let i = reference_base as usize;
                substitution_matrix[i][code] = read_base;
            }
        }

        Self(substitution_matrix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hit() {
        let mut frequencies = Frequencies::default();

        frequencies.hit(Base::A, Base::C);
        frequencies.hit(Base::G, Base::T);
        frequencies.hit(Base::T, Base::A);

        frequencies.hit(Base::C, Base::G);
        frequencies.hit(Base::C, Base::G);

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
    fn test_row() {
        let mut frequencies = Frequencies::default();
        frequencies.hit(Base::C, Base::A);
        assert_eq!(frequencies.row(Base::C), &[1, 0, 0, 0, 0]);
    }

    #[test]
    fn test_from_frequencies_for_substitution_matrix() {
        let frequencies = Frequencies([
            [0, 3, 8, 5, 0],
            [2, 0, 5, 5, 2],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        ]);

        let substitution_matrix = SubstitutionMatrix::from(frequencies);

        assert_eq!(
            substitution_matrix.0,
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
