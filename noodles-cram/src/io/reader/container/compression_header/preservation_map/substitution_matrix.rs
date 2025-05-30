use std::io;

use crate::container::compression_header::preservation_map::{
    SubstitutionMatrix, substitution_matrix::READ_BASES,
};

pub(super) fn read_substitution_matrix(src: &mut &[u8]) -> io::Result<SubstitutionMatrix> {
    let (buf, rest) = src
        .split_first_chunk()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(decode(*buf))
}

#[allow(clippy::identity_op)]
fn decode(src: [u8; 5]) -> SubstitutionMatrix {
    let mut substitution_matrix = SubstitutionMatrix::default();
    let matrix = &mut substitution_matrix.0;

    for ((substitutions, codes), read_bases) in matrix.iter_mut().zip(&src).zip(&READ_BASES) {
        substitutions[usize::from((codes >> 6) & 0x03)] = read_bases[0];
        substitutions[usize::from((codes >> 4) & 0x03)] = read_bases[1];
        substitutions[usize::from((codes >> 2) & 0x03)] = read_bases[2];
        substitutions[usize::from((codes >> 0) & 0x03)] = read_bases[3];
    }

    substitution_matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::container::compression_header::preservation_map::substitution_matrix::Base;

    #[test]
    fn test_read_substitution_matrix() -> io::Result<()> {
        // § 10.6.4 "Mapped reads: Substitution Matrix Format" (2024-09-04)
        let src = [0x93, 0x1b, 0x6c, 0xb1, 0xc6];
        let actual = read_substitution_matrix(&mut &src[..])?;

        let expected = SubstitutionMatrix([
            [Base::T, Base::G, Base::C, Base::N],
            [Base::A, Base::G, Base::T, Base::N],
            [Base::N, Base::A, Base::C, Base::T],
            [Base::G, Base::N, Base::A, Base::C],
            [Base::C, Base::G, Base::T, Base::A],
        ]);

        assert_eq!(actual, expected);

        Ok(())
    }
}
