use std::io::{self, Write};

use crate::file_definition::Version;

pub(super) fn write_format_version<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    let buf = [version.major(), version.minor()];
    writer.write_all(&buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_format_version() -> io::Result<()> {
        let mut buf = Vec::new();
        let version = Version::new(3, 1);
        write_format_version(&mut buf, version)?;
        assert_eq!(buf, [3, 1]);
        Ok(())
    }
}
