use std::io::{self, Write};

use crate::header::record::Kind;

pub(super) fn write_kind<W>(writer: &mut W, kind: Kind) -> io::Result<()>
where
    W: Write,
{
    let buf = match kind {
        Kind::Header => b"HD",
        Kind::ReferenceSequence => b"SQ",
        Kind::ReadGroup => b"RG",
        Kind::Program => b"PG",
        Kind::Comment => b"CO",
    };

    writer.write_all(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_kind() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, kind: Kind, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_kind(buf, kind)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Kind::Header, b"HD")?;
        t(&mut buf, Kind::ReferenceSequence, b"SQ")?;
        t(&mut buf, Kind::ReadGroup, b"RG")?;
        t(&mut buf, Kind::Program, b"PG")?;
        t(&mut buf, Kind::Comment, b"CO")?;

        Ok(())
    }
}
