use std::io::{self, Read};

use super::read_u8;

pub fn read_uint7<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut n = 0;

    loop {
        let b = read_u8(reader).map(u32::from)?;

        n <<= 7;
        n |= b & 0x7f;

        if b & 0x80 == 0 {
            break;
        }
    }

    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_uint7() -> io::Result<()> {
        fn t(mut data: &[u8], expected: u32) -> io::Result<()> {
            assert_eq!(read_uint7(&mut data)?, expected);
            Ok(())
        }

        // Examples from <https://en.wikipedia.org/wiki/Variable-length_quantity#Examples>.
        t(&[0x00], 0)?;
        t(&[0x7f], 127)?;
        t(&[0x81, 0x00], 128)?;
        t(&[0xc0, 0x00], 8192)?;
        t(&[0xff, 0x7f], 16383)?;
        t(&[0x81, 0x80, 0x00], 16384)?;
        t(&[0xff, 0xff, 0x7f], 2097151)?;
        t(&[0x81, 0x80, 0x80, 0x00], 2097152)?;
        t(&[0xc0, 0x80, 0x80, 0x00], 134217728)?;
        t(&[0xff, 0xff, 0xff, 0x7f], 268435455)?;

        Ok(())
    }
}
