use std::io::{self, Write};

use crate::io::writer::record::value::percent_encode;

pub(super) fn write_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    for t in percent_encode(s) {
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
