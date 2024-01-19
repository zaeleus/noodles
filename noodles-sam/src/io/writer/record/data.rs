mod field;

use std::io::{self, Write};

use self::field::write_field;
use crate::alignment::record::field::Data;

pub fn write_data<W, D>(writer: &mut W, data: D) -> io::Result<()>
where
    W: Write,
    D: Data,
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

        write_data(&mut buf, &data)?;

        assert_eq!(buf, b"\tNH:i:1\tCO:Z:noodles");

        Ok(())
    }
}
