use std::io::{self, Write};

use percent_encoding::{AsciiSet, CONTROLS, PercentEncode, utf8_percent_encode};

pub(super) fn write_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    for t in percent_encode(s) {
        writer.write_all(t.as_bytes())?;
    }

    Ok(())
}

// § 1.2 "Character encoding, non-printable characters and characters with special meaning"
// (2024-10-09).
const PERCENT_ENCODE_SET: &AsciiSet = &CONTROLS
    .add(b';')
    .add(b'=')
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
        let mut buf = Vec::new();

        buf.clear();
        write_string(&mut buf, "ndls")?;
        assert_eq!(buf, b"ndls");

        buf.clear();
        write_string(&mut buf, "n=d:ls;")?;
        assert_eq!(buf, b"n%3Dd:ls%3B");

        Ok(())
    }
}
