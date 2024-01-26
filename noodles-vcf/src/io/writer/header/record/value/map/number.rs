use std::io::{self, Write};

use crate::header::Number;

pub(super) fn write_number<W>(writer: &mut W, number: Number) -> io::Result<()>
where
    W: Write,
{
    const A: &[u8] = b"A";
    const R: &[u8] = b"R";
    const G: &[u8] = b"G";
    const UNKNOWN: &[u8] = b".";

    match number {
        Number::Count(n) => write!(writer, "{n}"),
        Number::A => writer.write_all(A),
        Number::R => writer.write_all(R),
        Number::G => writer.write_all(G),
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
        t(&mut buf, Number::A, b"A")?;
        t(&mut buf, Number::R, b"R")?;
        t(&mut buf, Number::G, b"G")?;
        t(&mut buf, Number::Unknown, b".")?;

        Ok(())
    }
}
