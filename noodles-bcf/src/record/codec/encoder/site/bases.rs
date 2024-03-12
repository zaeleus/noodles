use std::{
    io::{self, Write},
    str,
};

use noodles_vcf::variant::record::{AlternateBases, ReferenceBases};

use crate::record::codec::{encoder::value::write_value, Value};

pub(super) fn write_bases<W, R, A>(
    writer: &mut W,
    reference_bases: R,
    alternate_bases: A,
) -> io::Result<()>
where
    W: Write,
    R: ReferenceBases,
    A: AlternateBases,
{
    let raw_bases: Vec<_> = reference_bases.iter().collect::<io::Result<_>>()?;
    let r#ref =
        str::from_utf8(&raw_bases).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let ref_value = Some(Value::String(Some(r#ref)));
    write_value(writer, ref_value)?;

    if !alternate_bases.is_empty() {
        for result in alternate_bases.iter() {
            let alt = result?;
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
        use noodles_vcf::variant::record_buf::AlternateBases as AlternateBasesBuf;

        fn t(
            buf: &mut Vec<u8>,
            reference_bases: &str,
            alternate_bases: &AlternateBasesBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_bases(buf, reference_bases, alternate_bases)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "A", &AlternateBasesBuf::default(), &[0x17, b'A'])?;

        t(
            &mut buf,
            "A",
            &AlternateBasesBuf::from(vec![String::from("G")]),
            &[0x17, b'A', 0x17, b'G'],
        )?;

        t(
            &mut buf,
            "A",
            &AlternateBasesBuf::from(vec![String::from("G"), String::from("T")]),
            &[0x17, b'A', 0x17, b'G', 0x17, b'T'],
        )?;

        Ok(())
    }
}
