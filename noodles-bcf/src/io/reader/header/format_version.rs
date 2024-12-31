use std::io::{self, Read};

pub(super) fn read_format_version<R>(reader: &mut R) -> io::Result<(u8, u8)>
where
    R: Read,
{
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Ok((buf[0], buf[1]))
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
