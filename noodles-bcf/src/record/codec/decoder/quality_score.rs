use std::io;

use byteorder::{LittleEndian, ReadBytesExt};

pub(crate) fn read_qual(src: &mut &[u8]) -> io::Result<Option<f32>> {
    use crate::record::codec::value::Float;

    match src.read_f32::<LittleEndian>().map(Float::from)? {
        Float::Value(n) => Ok(Some(n)),
        Float::Missing => Ok(None),
        qual => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid qual: {qual:?}"),
        )),
    }
}
