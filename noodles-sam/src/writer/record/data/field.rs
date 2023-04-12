mod tag;
mod ty;
mod value;

use std::io::{self, Write};

use self::{tag::write_tag, ty::write_type, value::write_value};
use crate::record::data::field::{Tag, Value};

pub fn write_field<W>(writer: &mut W, tag: Tag, value: &Value) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: u8 = b':';

    write_tag(writer, tag)?;
    writer.write_all(&[DELIMITER])?;
    write_type(writer, value.ty())?;
    writer.write_all(&[DELIMITER])?;
    write_value(writer, value)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let (tag, value) = (Tag::AlignmentHitCount, Value::from(1));
        write_field(&mut buf, tag, &value)?;
        assert_eq!(buf, b"NH:i:1");
        Ok(())
    }
}
