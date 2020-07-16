mod base;

pub use self::base::Base;

use std::{convert::TryFrom, error, fmt};

type Substitutions = [[Base; 4]; 5];

#[derive(Debug, Default)]
pub struct SubstitutionMatrix {
    substitutions: Substitutions,
}

impl SubstitutionMatrix {
    pub fn get(&self, reference_base: Base, substitution_code: u8) -> Base {
        self.substitutions[reference_base as usize][substitution_code as usize]
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
}
