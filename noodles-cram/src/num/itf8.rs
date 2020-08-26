use std::io::{self, Read, Write};

use byteorder::{ReadBytesExt, WriteBytesExt};

use super::Itf8;

fn read_u8_as_i32<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_u8().map(i32::from)
}

pub fn read_itf8<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    let b0 = read_u8_as_i32(reader)?;

    let value = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        let b1 = read_u8_as_i32(reader)?;
        (b0 & 0x7f) << 8 | b1
    } else if b0 & 0x20 == 0 {
        let b1 = read_u8_as_i32(reader)?;
        let b2 = read_u8_as_i32(reader)?;
        (b0 & 0x3f) << 16 | b1 << 8 | b2
    } else if b0 & 0x10 == 0 {
        let b1 = read_u8_as_i32(reader)?;
        let b2 = read_u8_as_i32(reader)?;
        let b3 = read_u8_as_i32(reader)?;
        (b0 & 0x1f) << 24 | b1 << 16 | b2 << 8 | b3
    } else {
        let b1 = read_u8_as_i32(reader)?;
        let b2 = read_u8_as_i32(reader)?;
        let b3 = read_u8_as_i32(reader)?;
        let b4 = read_u8_as_i32(reader)?;
        (b0 & 0x0f) << 28 | b1 << 20 | b2 << 12 | b3 << 4 | b4 & 0x0f
    };

    Ok(value)
}

pub fn write_itf8<W>(writer: &mut W, value: Itf8) -> io::Result<()>
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

pub fn size_of(value: Itf8) -> usize {
    if value >> (8 - 1) == 0 {
        1
    } else if value >> (16 - 2) == 0 {
        2
    } else if value >> (24 - 3) == 0 {
        3
    } else if value >> (32 - 4) == 0 {
        4
    } else {
        5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_itf8() -> io::Result<()> {
        let data = [0x00];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 0);

        let data = [0x87, 0x55];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 1877);

        let data = [0xc7, 0x55, 0x99];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 480665);

        let data = [0xe7, 0x55, 0x99, 0x66];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 123050342);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x02];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 1968805474);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x12];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 1968805474);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x22];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 1968805474);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x42];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 1968805474);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x82];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, 1968805474);

        let data = [0xff, 0xff, 0xff, 0xff, 0x0f];
        let mut reader = &data[..];
        assert_eq!(read_itf8(&mut reader)?, -1);

        Ok(())
    }

    #[test]
    fn test_write_itf8() -> io::Result<()> {
        let mut buf = Vec::new();
        write_itf8(&mut buf, 0)?;
        assert_eq!(buf, [0x00]);

        let mut buf = Vec::new();
        write_itf8(&mut buf, 1877)?;
        assert_eq!(buf, [0x87, 0x55]);

        let mut buf = Vec::new();
        write_itf8(&mut buf, 480665)?;
        assert_eq!(buf, [0xc7, 0x55, 0x99]);

        let mut buf = Vec::new();
        write_itf8(&mut buf, 123050342)?;
        assert_eq!(buf, [0xe7, 0x55, 0x99, 0x66]);

        let mut buf = Vec::new();
        write_itf8(&mut buf, 1968805474)?;
        assert_eq!(buf, [0xf7, 0x55, 0x99, 0x66, 0x02]);

        let mut buf = Vec::new();
        write_itf8(&mut buf, -1)?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff, 0x0f]);

        Ok(())
    }

    #[test]
    fn test_size_of() {
        assert_eq!(size_of(0), 1);
        assert_eq!(size_of(1877), 2);
        assert_eq!(size_of(480665), 3);
        assert_eq!(size_of(123050342), 4);
        assert_eq!(size_of(1968805474), 5);
        assert_eq!(size_of(-1), 5);
    }
}
