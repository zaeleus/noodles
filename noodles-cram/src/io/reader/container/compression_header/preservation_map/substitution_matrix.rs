use std::io;

use crate::{
    container::compression_header::preservation_map::SubstitutionMatrix,
    io::reader::split_first_chunk,
};

pub(super) fn read_substitution_matrix(src: &mut &[u8]) -> io::Result<SubstitutionMatrix> {
    let (buf, rest) =
        split_first_chunk(src).ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    SubstitutionMatrix::try_from(*buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
