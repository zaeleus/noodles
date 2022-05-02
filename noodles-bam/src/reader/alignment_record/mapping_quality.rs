use std::{io, mem};

use bytes::Buf;
use noodles_sam::{self as sam, alignment::record::MappingQuality};

pub fn get_mapping_quality<B>(src: &mut B) -> io::Result<Option<MappingQuality>>
where
    B: Buf,
{
    use sam::alignment::record::mapping_quality::MISSING;

    if src.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_u8() {
        MISSING => Ok(None),
        n => Ok(MappingQuality::new(n)),
    }
}
