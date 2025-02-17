mod base;

use std::{error, fmt};

pub use self::base::Base;

type Substitutions = [[Base; 4]; 5];

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SubstitutionMatrix {
    pub(crate) substitutions: Substitutions,
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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromByteArrayError([u8; 5]);

impl fmt::Display for TryFromByteArrayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid substitution matrix: {:#x?}", self.0)
    }
}

impl error::Error for TryFromByteArrayError {}

impl TryFrom<[u8; 5]> for SubstitutionMatrix {
    type Error = TryFromByteArrayError;

    fn try_from(b: [u8; 5]) -> Result<Self, Self::Error> {
        let mut matrix = Self::default();

        set_substitutions(
            Base::A,
            b[0],
            [Base::C, Base::G, Base::T, Base::N],
            &mut matrix.substitutions,
        );

        set_substitutions(
            Base::C,
            b[1],
            [Base::A, Base::G, Base::T, Base::N],
            &mut matrix.substitutions,
        );

        set_substitutions(
            Base::G,
            b[2],
            [Base::A, Base::C, Base::T, Base::N],
            &mut matrix.substitutions,
        );

        set_substitutions(
            Base::T,
            b[3],
            [Base::A, Base::C, Base::G, Base::N],
            &mut matrix.substitutions,
        );

        set_substitutions(
            Base::N,
            b[4],
            [Base::A, Base::C, Base::G, Base::T],
            &mut matrix.substitutions,
        );

        Ok(matrix)
    }
}

fn set_substitutions(
    reference_base: Base,
    codes: u8,
    read_bases: [Base; 4],
    substitutions: &mut Substitutions,
) {
    let s = &mut substitutions[reference_base as usize];
    s[((codes >> 6) & 0x03) as usize] = read_bases[0];
    s[((codes >> 4) & 0x03) as usize] = read_bases[1];
    s[((codes >> 2) & 0x03) as usize] = read_bases[2];
    s[((codes) & 0x03) as usize] = read_bases[3];
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
    fn test_try_from_u8_slice() -> Result<(), TryFromByteArrayError> {
        let codes = [0x93, 0x1b, 0x6c, 0xb1, 0xc6];
        let matrix = SubstitutionMatrix::try_from(codes)?;

        let actual = &matrix.substitutions;
        let expected = &[
            [Base::T, Base::G, Base::C, Base::N],
            [Base::A, Base::G, Base::T, Base::N],
            [Base::N, Base::A, Base::C, Base::T],
            [Base::G, Base::N, Base::A, Base::C],
            [Base::C, Base::G, Base::T, Base::A],
        ];

        assert_eq!(actual, expected);

        Ok(())
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
