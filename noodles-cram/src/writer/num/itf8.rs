use std::io::{self, Write};

use byteorder::WriteBytesExt;

pub fn write_itf8<W>(writer: &mut W, value: i32) -> io::Result<()>
where
    W: Write,
{
    if value >> (8 - 1) == 0 {
        writer.write_u8(value as u8)?;
    } else if value >> (16 - 2) == 0 {
        writer.write_u8((value >> 8 | 0x80) as u8)?;
        writer.write_u8(value as u8)?;
    } else if value >> (24 - 3) == 0 {
        writer.write_u8((value >> 16 | 0xc0) as u8)?;
        writer.write_u8((value >> 8) as u8)?;
        writer.write_u8(value as u8)?;
    } else if value >> (32 - 4) == 0 {
        writer.write_u8((value >> 24 | 0xe0) as u8)?;
        writer.write_u8((value >> 16) as u8)?;
        writer.write_u8((value >> 8) as u8)?;
        writer.write_u8(value as u8)?;
    } else {
        writer.write_u8((value >> 28 | 0xf0) as u8)?;
        writer.write_u8((value >> 20) as u8)?;
        writer.write_u8((value >> 12) as u8)?;
        writer.write_u8((value >> 4) as u8)?;
        // § 2.3 Writing bytes to a byte stream: "only [the] 4 lower bits [are] used in the last
        // byte 5."
        writer.write_u8((value & 0x0f) as u8)?;
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
