use std::io::{self, Write};

use crate::header::record::value::map::format::Number;

pub(super) fn write_number<W>(writer: &mut W, number: Number) -> io::Result<()>
where
    W: Write,
{
    const ALTERNATE_BASES: &[u8] = b"A";
    const REFERENCE_ALTERNATE_BASES: &[u8] = b"R";
    const SAMPLES: &[u8] = b"G";
    const UNKNOWN: &[u8] = b".";

    match number {
        Number::Count(n) => write!(writer, "{n}"),
        Number::AlternateBases => writer.write_all(ALTERNATE_BASES),
        Number::ReferenceAlternateBases => writer.write_all(REFERENCE_ALTERNATE_BASES),
        Number::Samples => writer.write_all(SAMPLES),
        Number::Unknown => writer.write_all(UNKNOWN),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_number() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, number: Number, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_number(buf, number)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Number::Count(1), b"1")?;
        t(&mut buf, Number::AlternateBases, b"A")?;
        t(&mut buf, Number::ReferenceAlternateBases, b"R")?;
        t(&mut buf, Number::Samples, b"G")?;
        t(&mut buf, Number::Unknown, b".")?;

        Ok(())
    }
}
