//! BAM record data field component readers.

mod tag;
mod value;

pub use self::value::read_value;

use std::io::{self, BufRead};

use crate::record::data::Field;

pub(crate) fn read_field<R>(reader: &mut R) -> io::Result<Option<Field>>
where
    R: BufRead,
{
    use self::tag::read_tag;

    let tag = match read_tag(reader) {
        Ok(t) => t,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    };

    let ty = value::read_type(reader)?;
    let value = read_value(reader, ty)?;

    Ok(Some(Field::new(tag, value)))
}
