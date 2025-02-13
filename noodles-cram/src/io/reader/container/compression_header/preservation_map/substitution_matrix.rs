use std::io;

use crate::{
    container::compression_header::preservation_map::SubstitutionMatrix,
    io::reader::split_at_checked,
};

const SIZE: usize = 5;

pub(super) fn read_substitution_matrix(src: &mut &[u8]) -> io::Result<SubstitutionMatrix> {
    let (buf, rest) =
        split_at_checked(src, SIZE).ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    // SAFETY: `buf.len() == 5`.
    let encoded_matrix: [u8; SIZE] = buf.try_into().unwrap();

    SubstitutionMatrix::try_from(encoded_matrix)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
