//! BAM record data field component readers.

mod tag;
mod value;

pub use self::value::get_value;

use std::io;

use bytes::Buf;

use crate::record::data::Field;

pub(crate) fn get_field<B>(src: &mut B) -> io::Result<Option<Field>>
where
    B: Buf,
{
    use self::tag::get_tag;

    let tag = match get_tag(src) {
        Ok(t) => t,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    };

    let ty = value::get_type(src)?;
    let value = get_value(src, ty)?;

    Ok(Some(Field::new(tag, value)))
}
