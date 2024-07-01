use std::io::{self, Write};

use crate::header::record::value::map::format::Number;

pub(super) fn write_number<W>(writer: &mut W, number: Number) -> io::Result<()>
where
    W: Write,
{
    match number {
        Number::Count(n) => write!(writer, "{n}"),
        Number::AlternateBases => writer.write_all(b"A"),
        Number::ReferenceAlternateBases => writer.write_all(b"R"),
        Number::Samples => writer.write_all(b"G"),
        Number::LocalAlternateBases => writer.write_all(b"LA"),
        Number::LocalReferenceAlternateBases => writer.write_all(b"LR"),
        Number::LocalSamples => writer.write_all(b"LG"),
        Number::Ploidy => writer.write_all(b"P"),
        Number::BaseModifications => writer.write_all(b"M"),
        Number::Unknown => writer.write_all(b"."),
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
        t(&mut buf, Number::LocalAlternateBases, b"LA")?;
        t(&mut buf, Number::LocalReferenceAlternateBases, b"LR")?;
        t(&mut buf, Number::LocalSamples, b"LG")?;
        t(&mut buf, Number::Ploidy, b"P")?;
        t(&mut buf, Number::BaseModifications, b"M")?;
        t(&mut buf, Number::Unknown, b".")?;

        Ok(())
    }
}
