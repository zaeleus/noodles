use std::io::{self, Write};

use super::write_missing;
use crate::record_buf::Strand;

pub(super) fn write_strand<W>(writer: &mut W, strand: Option<Strand>) -> io::Result<()>
where
    W: Write,
{
    const FORWARD: &[u8] = b"+";
    const REVERSE: &[u8] = b"-";

    match strand {
        Some(Strand::Forward) => writer.write_all(FORWARD),
        Some(Strand::Reverse) => writer.write_all(REVERSE),
        None => write_missing(writer),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_strand() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, strand: Option<Strand>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_strand(buf, strand)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Some(Strand::Forward), b"+")?;
        t(&mut buf, Some(Strand::Reverse), b"-")?;
        t(&mut buf, None, b".")?;

        Ok(())
    }
}
