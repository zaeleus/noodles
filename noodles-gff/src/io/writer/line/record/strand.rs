use std::io::{self, Write};

use super::write_missing;
use crate::record_buf::Strand;

const FORWARD: &[u8] = b"+";
const REVERSE: &[u8] = b"-";
const UNKNOWN: &[u8] = b"?";

pub(super) fn write_strand<W>(writer: &mut W, strand: Strand) -> io::Result<()>
where
    W: Write,
{
    match strand {
        Strand::None => write_missing(writer),
        Strand::Forward => writer.write_all(FORWARD),
        Strand::Reverse => writer.write_all(REVERSE),
        Strand::Unknown => writer.write_all(UNKNOWN),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_strand() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, strand: Strand, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_strand(buf, strand)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Strand::None, b".")?;
        t(&mut buf, Strand::Forward, b"+")?;
        t(&mut buf, Strand::Reverse, b"-")?;
        t(&mut buf, Strand::Unknown, b"?")?;

        Ok(())
    }
}
