use std::io::{self, Write};

use percent_encoding::{utf8_percent_encode, AsciiSet, CONTROLS};

// § 1.2 "Character encoding, non-printable characters and characters with special meaning" (2023-08-23)
const PERCENT_ENCODE_SET: &AsciiSet = &CONTROLS
    .add(b':')
    .add(b';')
    .add(b'=')
    .add(b'%')
    .add(b',')
    .add(b'\r')
    .add(b'\n')
    .add(b'\t');

pub(super) fn write_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    for t in utf8_percent_encode(s, PERCENT_ENCODE_SET) {
        writer.write_all(t.as_bytes())?;
    }

    Ok(())
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
        write_string(&mut buf, "noodles=vcf;")?;
        assert_eq!(buf, b"noodles%3Dvcf%3B");

        Ok(())
    }
}
