use std::io;

use bytes::Buf;

use crate::container::compression_header::preservation_map::SubstitutionMatrix;

pub(super) fn get_substitution_matrix<B>(src: &mut B) -> io::Result<SubstitutionMatrix>
where
    B: Buf,
{
    let mut buf = [0; 5];

    if src.remaining() < buf.len() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    src.copy_to_slice(&mut buf);

    SubstitutionMatrix::try_from(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
