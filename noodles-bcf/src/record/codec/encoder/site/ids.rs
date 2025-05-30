use std::io::{self, Write};

use noodles_vcf::variant::record::Ids;

use crate::record::codec::{Value, encoder::value::write_value};

pub(super) fn write_ids<W, I>(writer: &mut W, ids: I) -> io::Result<()>
where
    W: Write,
    I: Ids,
{
    const DELIMITER: &str = ";";

    if ids.is_empty() {
        let value = Some(Value::String(None));
        write_value(writer, value)
    } else {
        let s = ids.iter().collect::<Vec<_>>().join(DELIMITER);
        let value = Some(Value::String(Some(&s)));
        write_value(writer, value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_ids() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_vcf::variant::record_buf::Ids as IdsBuf;

        fn t(buf: &mut Vec<u8>, ids: &IdsBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_ids(buf, ids)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &IdsBuf::default(), &[0x07])?;
        t(
            &mut buf,
            &[String::from("nd0")].into_iter().collect(),
            &[0x37, b'n', b'd', b'0'],
        )?;
        t(
            &mut buf,
            &[String::from("nd0"), String::from("nd1")]
                .into_iter()
                .collect(),
            &[0x77, b'n', b'd', b'0', b';', b'n', b'd', b'1'],
        )?;

        Ok(())
    }
}
