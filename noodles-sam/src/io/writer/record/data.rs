mod field;
mod typed_map_encoded;

use std::io::{self, Write};

use self::{field::write_field, typed_map_encoded::write_typed_map_encoded_data};
use crate::alignment::record::{Data, DataRef};

const SEPARATOR: u8 = b'\t';

pub(super) fn write_data<W>(writer: &mut W, data: DataRef<'_>) -> io::Result<()>
where
    W: Write,
{
    match data {
        DataRef::TypedMapEncoded(src) => write_typed_map_encoded_data(writer, src)?,
        DataRef::Data(data) => write_generic_data(writer, data)?,
    }

    Ok(())
}

pub(super) fn write_generic_data<'r, W, D>(writer: &mut W, data: D) -> io::Result<()>
where
    W: Write,
    D: Data<'r>,
{
    for result in data.iter() {
        let (tag, value) = result?;

        writer.write_all(&[SEPARATOR])?;
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

        let d = DataRef::Data(Box::new(&data));
        write_data(&mut buf, d)?;

        assert_eq!(buf, b"\tNH:i:1\tCO:Z:noodles");

        Ok(())
    }
}
