use std::io::{self, Write};

use percent_encoding::{
    AsciiSet, CONTROLS, PercentEncode, percent_encode_byte, utf8_percent_encode,
};

const MISSING: &[u8; 1] = b".";

pub(super) fn write_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    if s.as_bytes() == MISSING {
        let t = percent_encode_byte(MISSING[0]);
        writer.write_all(t.as_bytes())?;
    } else {
        for t in percent_encode(s) {
            writer.write_all(t.as_bytes())?;
        }
    }

    Ok(())
}

// § 1.2 "Character encoding, non-printable characters and characters with special meaning"
// (2024-10-09).
const PERCENT_ENCODE_SET: &AsciiSet = &CONTROLS
    .add(b':')
    .add(b'%')
    .add(b',')
    .add(b'\r')
    .add(b'\n')
    .add(b'\t');

pub(super) fn percent_encode(s: &str) -> PercentEncode<'_> {
    utf8_percent_encode(s, PERCENT_ENCODE_SET)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_string() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, s: &str, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_string(buf, s)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "ndls", b"ndls")?;
        t(&mut buf, "n=d:ls.;", b"n=d%3Als.;")?;
        t(&mut buf, ".", b"%2E")?;

        Ok(())
    }
}
