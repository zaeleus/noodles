mod base;

pub use self::base::Base;

pub const READ_BASES: [[Base; 4]; 5] = [
    [Base::C, Base::G, Base::T, Base::N], // A
    [Base::A, Base::G, Base::T, Base::N], // C
    [Base::A, Base::C, Base::T, Base::N], // G
    [Base::A, Base::C, Base::G, Base::N], // T
    [Base::A, Base::C, Base::G, Base::T], // N
];

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SubstitutionMatrix {
    pub(crate) substitutions: [[Base; 4]; 5],
}

impl SubstitutionMatrix {
    pub fn get(&self, reference_base: Base, substitution_code: u8) -> Base {
        let i = match reference_base {
            Base::A => 0,
            Base::C => 1,
            Base::G => 2,
            Base::T => 3,
            Base::N => 4,
        };

        let j = usize::from(substitution_code & 0x03);

        self.substitutions[i][j]
    }

    pub fn find_code(&self, reference_base: Base, read_base: Base) -> u8 {
        for code in [0b00, 0b01, 0b10, 0b11] {
            if self.get(reference_base, code) == read_base {
                return code;
            }
        }

        unreachable!();
    }
}

impl Default for SubstitutionMatrix {
    fn default() -> Self {
        SubstitutionMatrix {
            substitutions: [
                [Base::C, Base::G, Base::T, Base::N],
                [Base::A, Base::G, Base::T, Base::N],
                [Base::A, Base::C, Base::T, Base::N],
                [Base::A, Base::C, Base::G, Base::N],
                [Base::A, Base::C, Base::G, Base::T],
            ],
        }
    }
}

impl From<SubstitutionMatrix> for [u8; 5] {
    fn from(substitution_matrix: SubstitutionMatrix) -> Self {
        <[u8; 5]>::from(&substitution_matrix)
    }
}

impl From<&SubstitutionMatrix> for [u8; 5] {
    fn from(substitution_matrix: &SubstitutionMatrix) -> Self {
        let mut encoded_substitution_matrix = [0; 5];

        for (read_bases, codes) in substitution_matrix
            .substitutions
            .iter()
            .zip(encoded_substitution_matrix.iter_mut())
        {
            let mut index_base_pairs: Vec<_> = read_bases.iter().enumerate().collect();
            index_base_pairs.sort_by_key(|p| p.1);

            for (i, _) in index_base_pairs {
                *codes <<= 2;
                *codes |= i as u8;
            }
        }

        encoded_substitution_matrix
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_code() {
        let matrix = SubstitutionMatrix {
            substitutions: [
                [Base::T, Base::G, Base::C, Base::N],
                [Base::A, Base::G, Base::T, Base::N],
                [Base::N, Base::A, Base::C, Base::T],
                [Base::G, Base::N, Base::A, Base::C],
                [Base::C, Base::G, Base::T, Base::A],
            ],
        };

        assert_eq!(matrix.find_code(Base::A, Base::T), 0b00);
        assert_eq!(matrix.find_code(Base::C, Base::G), 0b01);
        assert_eq!(matrix.find_code(Base::G, Base::C), 0b10);
        assert_eq!(matrix.find_code(Base::T, Base::C), 0b11);
    }

    #[test]
    fn test_from_substitution_matrix_for_5_byte_array() {
        let matrix = SubstitutionMatrix {
            substitutions: [
                [Base::T, Base::G, Base::C, Base::N],
                [Base::A, Base::G, Base::T, Base::N],
                [Base::N, Base::A, Base::C, Base::T],
                [Base::G, Base::N, Base::A, Base::C],
                [Base::C, Base::G, Base::T, Base::A],
            ],
        };

        assert_eq!(<[u8; 5]>::from(&matrix), [0x93, 0x1b, 0x6c, 0xb1, 0xc6]);
        assert_eq!(<[u8; 5]>::from(matrix), [0x93, 0x1b, 0x6c, 0xb1, 0xc6]);
    }
}
