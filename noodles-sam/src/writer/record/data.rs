mod field;

use std::io::{self, Write};

use self::field::write_field;
use crate::record::Data;

pub fn write_data<W>(writer: &mut W, data: &Data) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: u8 = b'\t';

    for (tag, value) in data.iter() {
        writer.write_all(&[DELIMITER])?;
        write_field(writer, tag, value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_data() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::data::field::{tag, Value};

        let mut buf = Vec::new();

        let data = [
            (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
            (tag::COMMENT, Value::try_from(String::from("noodles"))?),
        ]
        .into_iter()
        .collect();

        write_data(&mut buf, &data)?;

        assert_eq!(buf, b"\tNH:i:1\tCO:Z:noodles");

        Ok(())
    }
}
