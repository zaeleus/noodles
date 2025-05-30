use std::io::{self, Write};

use crate::header::record::value::map::header::{Version, tag};

pub(super) fn write_version_field<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    use crate::io::writer::header::record::{value::map::write_separator, write_delimiter};

    write_delimiter(writer)?;
    writer.write_all(tag::VERSION.as_ref())?;
    write_separator(writer)?;
    write_version(writer, version)?;

    Ok(())
}

fn write_version<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    use crate::io::writer::num;

    const DELIMITER: u8 = b'.';

    num::write_u32(writer, version.major())?;
    writer.write_all(&[DELIMITER])?;
    num::write_u32(writer, version.minor())?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_version_field() -> io::Result<()> {
        let mut buf = Vec::new();
        write_version_field(&mut buf, Version::new(1, 6))?;
        assert_eq!(buf, b"\tVN:1.6");
        Ok(())
    }
}
