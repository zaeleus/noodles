use std::io::{self, Write};

pub fn write_itf8<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    if n >> (8 - 1) == 0 {
        let buf = [n as u8];
        writer.write_all(&buf)?;
    } else if n >> (16 - 2) == 0 {
        let buf = [((n >> 8) | 0x80) as u8, n as u8];
        writer.write_all(&buf)?;
    } else if n >> (24 - 3) == 0 {
        let buf = [((n >> 16) | 0xc0) as u8, (n >> 8) as u8, n as u8];
        writer.write_all(&buf)?;
    } else if n >> (32 - 4) == 0 {
        let buf = [
            ((n >> 24) | 0xe0) as u8,
            (n >> 16) as u8,
            (n >> 8) as u8,
            n as u8,
        ];

        writer.write_all(&buf)?;
    } else {
        let buf = [
            ((n >> 28) | 0xf0) as u8,
            (n >> 20) as u8,
            (n >> 12) as u8,
            (n >> 4) as u8,
            // § 2.3.4 "Writing bytes to a byte stream: ITF-8 integer (itf8)" (2024-09-04):
            // "...only [the] 4 lower bits [are] used in the last byte 5."
            (n & 0x0f) as u8,
        ];

        writer.write_all(&buf)?;
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
