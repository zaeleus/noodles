mod tag;
mod value;

pub use self::{tag::read_tag, value::read_value};

use std::io::{self, BufRead};

use crate::record::data::Field;

pub fn read_field<R>(reader: &mut R) -> io::Result<Field>
where
    R: BufRead,
{
    let tag = read_tag(reader)?;
    let ty = value::read_type(reader)?;
    let value = read_value(reader, ty)?;
    Ok(Field::new(tag, value))
}
