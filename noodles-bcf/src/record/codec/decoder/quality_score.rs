use std::io;

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::record::QualityScore;

pub(crate) fn read_qual(src: &mut &[u8]) -> io::Result<Option<QualityScore>> {
    use crate::lazy::record::value::Float;

    match src.read_f32::<LittleEndian>().map(Float::from)? {
        Float::Value(value) => QualityScore::try_from(value)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)),
        Float::Missing => Ok(None),
        qual => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid qual: {qual:?}"),
        )),
    }
}
