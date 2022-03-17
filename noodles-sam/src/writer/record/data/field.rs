mod tag;
mod value;

use std::io::{self, Write};

use self::{tag::write_tag, value::write_value};
use crate::record::data::Field;

pub fn write_field<W>(writer: &mut W, field: &Field) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: u8 = b':';

    write_tag(writer, field.tag())?;
    writer.write_all(&[DELIMITER])?;
    value::write_type(writer, field.value().ty())?;
    writer.write_all(&[DELIMITER])?;
    write_value(writer, field.value())?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_field() -> io::Result<()> {
        use crate::record::data::field::{Tag, Value};

        let mut buf = Vec::new();
        let field = Field::new(Tag::AlignmentHitCount, Value::from(1));
        write_field(&mut buf, &field)?;
        assert_eq!(buf, b"NH:i:1");

        Ok(())
    }
}
