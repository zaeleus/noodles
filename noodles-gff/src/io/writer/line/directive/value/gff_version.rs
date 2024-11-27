use std::io::{self, Write};

use crate::directive_buf::value::GffVersion;

pub(crate) fn write_gff_version<W>(writer: &mut W, version: &GffVersion) -> io::Result<()>
where
    W: Write,
{
    write!(writer, "{}", version.major())?;

    if let Some(minor) = version.minor() {
        write_separator(writer)?;
        write!(writer, "{minor}")?;

        if let Some(patch) = version.patch() {
            write_separator(writer)?;
            write!(writer, "{patch}")?;
        }
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b'.';
    writer.write_all(&[SEPARATOR])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_gff_version() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, version: &GffVersion, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_gff_version(buf, version)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &"3".parse()?, b"3")?;
        t(&mut buf, &"3.1".parse()?, b"3.1")?;
        t(&mut buf, &"3.1.26".parse()?, b"3.1.26")?;

        Ok(())
    }
}
