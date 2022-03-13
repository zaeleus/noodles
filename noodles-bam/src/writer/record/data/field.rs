mod tag;
mod value;

use std::io::{self, Write};

use noodles_sam::record::data::Field;

use self::{tag::write_tag, value::write_value};

pub fn write_field<W>(writer: &mut W, field: &Field) -> io::Result<()>
where
    W: Write,
{
    write_tag(writer, field.tag())?;
    value::write_type(writer, field.value().ty())?;
    write_value(writer, field.value())?;
    Ok(())
}
