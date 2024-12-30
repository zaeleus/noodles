use std::io::{self, Read};

use byteorder::ReadBytesExt;

pub(crate) fn read_format_version<R>(reader: &mut R) -> io::Result<(u8, u8)>
where
    R: Read,
{
    let major_version = reader.read_u8()?;
    let minor_version = reader.read_u8()?;

    Ok((major_version, minor_version))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_format_version() -> io::Result<()> {
        let data = [0x02, 0x01];
        let mut reader = &data[..];
        assert_eq!(read_format_version(&mut reader)?, (2, 1));
        Ok(())
    }
}
