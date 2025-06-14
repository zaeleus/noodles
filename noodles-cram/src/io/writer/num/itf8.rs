use std::io::{self, Write};

pub fn write_itf8<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    if n >> (8 - 1) == 0 {
        let buf = [n as u8];
        writer.write_all(&buf)?;
    } else if n >> (16 - 2) == 0 {
        let m = (n as u16) | (0x80 << 8);
        let buf = m.to_be_bytes();
        writer.write_all(&buf)?;
    } else if n >> (24 - 3) == 0 {
        let m = (n as u32) | (0xc0 << 16);
        let buf = m.to_be_bytes();
        writer.write_all(&buf[1..])?;
    } else if n >> (32 - 4) == 0 {
        let m = (n as u32) | (0xe0 << 24);
        let buf = m.to_be_bytes();
        writer.write_all(&buf)?;
    } else {
        let m = n as u32 as u64;

        // § 2.3.4 "Writing bytes to a byte stream: ITF-8 integer (itf8)" (2024-09-04): "...only
        // [the] 4 lower bits [are] used in the last byte 5."
        let m = (0xf0 << 32) | ((m & 0xfffffff0) << 4) | (m & 0x0f);

        let buf = m.to_be_bytes();
        writer.write_all(&buf[3..])?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_itf8() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, value: i32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_itf8(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &[0x00])?;
        t(&mut buf, 1877, &[0x87, 0x55])?;
        t(&mut buf, 480665, &[0xc7, 0x55, 0x99])?;
        t(&mut buf, 123050342, &[0xe7, 0x55, 0x99, 0x66])?;
        t(&mut buf, 1968805474, &[0xf7, 0x55, 0x99, 0x66, 0x02])?;
        t(&mut buf, -1, &[0xff, 0xff, 0xff, 0xff, 0x0f])?;

        Ok(())
    }
}
