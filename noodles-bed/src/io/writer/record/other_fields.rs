mod value;

use std::io::{self, Write};

use self::value::write_value;
use super::write_separator;
use crate::feature::record::OtherFields;

pub(super) fn write_other_fields<W, F>(writer: &mut W, other_fields: &F) -> io::Result<()>
where
    W: Write,
    F: OtherFields + ?Sized,
{
    for value in other_fields.iter() {
        write_separator(writer)?;
        write_value(writer, value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::record_buf::{other_fields::Value, OtherFields as OtherFieldsBuf};

    #[test]
    fn test_write_other_fields() -> io::Result<()> {
        let mut buf = Vec::new();

        let mut other_fields = OtherFieldsBuf::default();
        buf.clear();
        write_other_fields(&mut buf, &other_fields)?;
        assert!(buf.is_empty());

        other_fields.as_mut().push(Value::UInt64(8));
        buf.clear();
        write_other_fields(&mut buf, &other_fields)?;
        assert_eq!(buf, b"\t8");

        other_fields.as_mut().push(Value::UInt64(13));
        buf.clear();
        write_other_fields(&mut buf, &other_fields)?;
        assert_eq!(buf, b"\t8\t13");

        Ok(())
    }
}
