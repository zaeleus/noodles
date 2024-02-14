use std::io::{self, Write};

use noodles_vcf as vcf;

use crate::record::codec::{encoder::value::write_value, Value};

pub(super) fn write_ids<W>(writer: &mut W, ids: &vcf::variant::record_buf::Ids) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &str = ";";

    if ids.as_ref().is_empty() {
        let value = Some(Value::String(None));
        write_value(writer, value)
    } else {
        let s = ids
            .as_ref()
            .iter()
            .map(|id| id.as_ref())
            .collect::<Vec<_>>()
            .join(DELIMITER);

        let value = Some(Value::String(Some(&s)));

        write_value(writer, value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_ids() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::variant::record_buf::Ids;

        fn t(buf: &mut Vec<u8>, ids: &Ids, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_ids(buf, ids)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Ids::default(), &[0x07])?;
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
