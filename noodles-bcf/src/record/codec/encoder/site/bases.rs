use std::io::{self, Write};

use noodles_vcf as vcf;

use crate::record::codec::{encoder::value::write_value, Value};

pub(super) fn write_bases<W>(
    writer: &mut W,
    reference_bases: &str,
    alternate_bases: &vcf::variant::record_buf::AlternateBases,
) -> io::Result<()>
where
    W: Write,
{
    let r#ref = reference_bases;
    let ref_value = Some(Value::String(Some(r#ref)));
    write_value(writer, ref_value)?;

    if !alternate_bases.is_empty() {
        for alt in alternate_bases.as_ref().iter() {
            let alt_value = Some(Value::String(Some(alt)));
            write_value(writer, alt_value)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_bases() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::variant::record_buf::AlternateBases;

        fn t(
            buf: &mut Vec<u8>,
            reference_bases: &str,
            alternate_bases: &AlternateBases,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_bases(buf, reference_bases, alternate_bases)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "A", &AlternateBases::default(), &[0x17, b'A'])?;

        t(
            &mut buf,
            "A",
            &AlternateBases::from(vec![String::from("G")]),
            &[0x17, b'A', 0x17, b'G'],
        )?;

        t(
            &mut buf,
            "A",
            &AlternateBases::from(vec![String::from("G"), String::from("T")]),
            &[0x17, b'A', 0x17, b'G', 0x17, b'T'],
        )?;

        Ok(())
    }
}
