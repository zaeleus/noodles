use std::io::{self, Write};

pub fn write_uint7<W>(writer: &mut W, mut n: u32) -> io::Result<()>
where
    W: Write,
{
    let mut buf = [0; 5];

    let mut i = buf.len() - 1;
    let b = (n & 0x7f) as u8;
    buf[i] = b;
    n >>= 7;

    while n > 0 {
        i -= 1;
        let b = (n & 0x7f) as u8;
        buf[i] = b | 0x80;
        n >>= 7;
    }

    writer.write_all(&buf[i..])?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_uint7() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: u32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_uint7(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        // Examples from <https://en.wikipedia.org/wiki/Variable-length_quantity#Examples>.
        t(&mut buf, 0, &[0x00])?;
        t(&mut buf, 127, &[0x7f])?;
        t(&mut buf, 128, &[0x81, 0x00])?;
        t(&mut buf, 8192, &[0xc0, 0x00])?;
        t(&mut buf, 16383, &[0xff, 0x7f])?;
        t(&mut buf, 16384, &[0x81, 0x80, 0x00])?;
        t(&mut buf, 2097151, &[0xff, 0xff, 0x7f])?;
        t(&mut buf, 2097152, &[0x81, 0x80, 0x80, 0x00])?;
        t(&mut buf, 134217728, &[0xc0, 0x80, 0x80, 0x00])?;
        t(&mut buf, 268435455, &[0xff, 0xff, 0xff, 0x7f])?;

        Ok(())
    }
}
