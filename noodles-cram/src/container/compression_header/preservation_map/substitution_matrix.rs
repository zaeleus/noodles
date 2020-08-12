mod base;
mod builder;
mod histogram;

pub use self::{base::Base, builder::Builder};

use std::{convert::TryFrom, error, fmt};

use self::histogram::Histogram;

type Substitutions = [[Base; 4]; 5];

#[derive(Debug)]
pub struct SubstitutionMatrix {
    substitutions: Substitutions,
}

impl SubstitutionMatrix {
    pub fn builder(reference_sequence: &[u8]) -> Builder<'_> {
        Builder::new(reference_sequence)
    }

    pub fn get(&self, reference_base: Base, substitution_code: u8) -> Base {
        self.substitutions[reference_base as usize][substitution_code as usize]
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

impl From<Histogram> for SubstitutionMatrix {
    fn from(histogram: Histogram) -> Self {
        let mut matrix = Self::default();
        let bases = [Base::A, Base::C, Base::G, Base::T, Base::N];

        for &reference_base in bases.iter() {
            let mut base_frequency_pairs: Vec<_> = bases
                .iter()
                .zip(histogram.get_bins(reference_base).iter())
                .filter(|(&read_base, _)| read_base != reference_base)
                .collect();

            // FIXME: u64 => i64 cast
            base_frequency_pairs.sort_by_key(|(_, &frequency)| -(frequency as i64));

            for (code, (&read_base, _)) in base_frequency_pairs.iter().enumerate() {
                let i = reference_base as usize;
                matrix.substitutions[i][code] = read_base;
            }
        }

        matrix
    }
}

#[derive(Debug)]
pub struct TryFromByteSliceError(Vec<u8>);

impl fmt::Display for TryFromByteSliceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid substitution matrix: {:#x?}", self.0)
    }
}

impl error::Error for TryFromByteSliceError {}

impl TryFrom<&[u8]> for SubstitutionMatrix {
    type Error = TryFromByteSliceError;

    fn try_from(bytes: &[u8]) -> Result<Self, Self::Error> {
        let mut matrix = Self::default();

        set_substitutions(
            Base::A,
            bytes
                .get(0)
                .copied()
                .ok_or_else(|| TryFromByteSliceError(bytes.to_vec()))?,
            [Base::C, Base::G, Base::T, Base::N],
            &mut matrix.substitutions,
        );

        set_substitutions(
            Base::C,
            bytes
                .get(1)
                .copied()
                .ok_or_else(|| TryFromByteSliceError(bytes.to_vec()))?,
            [Base::A, Base::G, Base::T, Base::N],
            &mut matrix.substitutions,
        );

        set_substitutions(
            Base::G,
            bytes
                .get(2)
                .copied()
                .ok_or_else(|| TryFromByteSliceError(bytes.to_vec()))?,
            [Base::A, Base::C, Base::T, Base::N],
            &mut matrix.substitutions,
        );

        set_substitutions(
            Base::T,
            bytes
                .get(3)
                .copied()
                .ok_or_else(|| TryFromByteSliceError(bytes.to_vec()))?,
            [Base::A, Base::C, Base::G, Base::N],
            &mut matrix.substitutions,
        );

        set_substitutions(
            Base::N,
            bytes
                .get(4)
                .copied()
                .ok_or_else(|| TryFromByteSliceError(bytes.to_vec()))?,
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
    fn test_try_from_u8_slice() -> Result<(), TryFromByteSliceError> {
        let codes = [0x93, 0x1b, 0x6c, 0xb1, 0xc6];
        let matrix = SubstitutionMatrix::try_from(&codes[..])?;

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

        assert_eq!(<[u8; 5]>::from(matrix), [0x93, 0x1b, 0x6c, 0xb1, 0xc6]);
    }

    #[test]
    fn test_from_histogram() {
        let histogram = Histogram::new([
            [0, 3, 8, 5, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        ]);

        let substitution_matrix = SubstitutionMatrix::from(histogram);

        let expected = [
            [Base::G, Base::T, Base::C, Base::N],
            [Base::A, Base::G, Base::T, Base::N],
            [Base::A, Base::C, Base::T, Base::N],
            [Base::A, Base::C, Base::G, Base::N],
            [Base::A, Base::C, Base::G, Base::T],
        ];

        assert_eq!(substitution_matrix.substitutions, expected);
    }
}
