//! BAM record data field component readers.

mod tag;
mod value;

pub use self::value::read_value;

use std::io::{self, BufRead};

use crate::record::data::Field;

pub(crate) fn read_field<R>(reader: &mut R) -> io::Result<Field>
where
    R: BufRead,
{
    use self::tag::read_tag;

    let tag = read_tag(reader)?;
    let ty = value::read_type(reader)?;
    let value = read_value(reader, ty)?;
    Ok(Field::new(tag, value))
}
