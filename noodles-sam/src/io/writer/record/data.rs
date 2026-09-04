mod field;
mod field_encoded;

use std::io::{self, Write};

use self::{field::write_field, field_encoded::write_field_encoded_data};
use crate::alignment::record::{Data, DataRef};

pub(super) fn write_data<W>(writer: &mut W, data: DataRef<'_>) -> io::Result<()>
where
    W: Write,
{
    match data {
        DataRef::FieldEncoded(src) => write_field_encoded_data(writer, src),
        DataRef::Data(d) => write_generic_data(writer, d),
    }
}

fn write_generic_data<'r, W, D>(writer: &mut W, data: D) -> io::Result<()>
where
    W: Write,
    D: Data<'r>,
{
    const DELIMITER: u8 = b'\t';

    for result in data.iter() {
        let (tag, value) = result?;

        writer.write_all(&[DELIMITER])?;
        write_field(writer, tag, &value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::record_buf::Data as DataBuf;

    #[test]
    fn test_write_data() -> io::Result<()> {
        use crate::alignment::{record::data::field::Tag, record_buf::data::field::Value};

        let mut buf = Vec::new();

        let data: DataBuf = [
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
            (Tag::COMMENT, Value::from("noodles")),
        ]
        .into_iter()
        .collect();

        write_data(&mut buf, DataRef::Data(Box::new(&data)))?;

        assert_eq!(buf, b"\tNH:i:1\tCO:Z:noodles");

        Ok(())
    }
}
