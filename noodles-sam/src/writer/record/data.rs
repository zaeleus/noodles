mod field;

use std::io::{self, Write};

use self::field::write_field;
use crate::record::Data;

pub fn write_data<W>(writer: &mut W, data: &Data) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: u8 = b'\t';

    for field in data.values() {
        writer.write_all(&[DELIMITER])?;
        write_field(writer, field)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_data() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::data::{
            field::{Tag, Value},
            Field,
        };

        let mut buf = Vec::new();

        let data = Data::try_from(vec![
            Field::new(Tag::AlignmentHitCount, Value::from(1)),
            Field::new(Tag::Comment, Value::try_from(String::from("noodles"))?),
        ])?;

        write_data(&mut buf, &data)?;

        assert_eq!(buf, b"\tNH:i:1\tCO:Z:noodles");

        Ok(())
    }
}
