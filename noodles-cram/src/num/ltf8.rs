use std::io::{self, Read, Write};

use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};

use super::Ltf8;

fn read_u8_as_i64<R>(reader: &mut R) -> io::Result<i64>
where
    R: Read,
{
    reader.read_u8().map(|b| b as i64)
}

pub fn read_ltf8<R>(reader: &mut R) -> io::Result<i64>
where
    R: Read,
{
    let b0 = read_u8_as_i64(reader)?;

    let value = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        let b1 = read_u8_as_i64(reader)?;
        (b0 & 0x7f) << 8 | b1
    } else if b0 & 0x20 == 0 {
        let b1 = read_u8_as_i64(reader)?;
        let b2 = read_u8_as_i64(reader)?;
        (b0 & 0x3f) << 16 | b1 << 8 | b2
    } else if b0 & 0x10 == 0 {
        let b1 = read_u8_as_i64(reader)?;
        let b2 = read_u8_as_i64(reader)?;
        let b3 = read_u8_as_i64(reader)?;
        (b0 & 0x1f) << 24 | b1 << 16 | b2 << 8 | b3
    } else if b0 & 0x08 == 0 {
        let b1 = read_u8_as_i64(reader)?;
        let b2 = read_u8_as_i64(reader)?;
        let b3 = read_u8_as_i64(reader)?;
        let b4 = read_u8_as_i64(reader)?;
        (b0 & 0x0f) << 32 | b1 << 24 | b2 << 16 | b3 << 8 | b4
    } else if b0 & 0x04 == 0 {
        let b1 = read_u8_as_i64(reader)?;
        let b2 = read_u8_as_i64(reader)?;
        let b3 = read_u8_as_i64(reader)?;
        let b4 = read_u8_as_i64(reader)?;
        let b5 = read_u8_as_i64(reader)?;
        (b0 & 0x07) << 40 | b1 << 32 | b2 << 24 | b3 << 16 | b4 << 8 | b5
    } else if b0 & 0x02 == 0 {
        let b1 = read_u8_as_i64(reader)?;
        let b2 = read_u8_as_i64(reader)?;
        let b3 = read_u8_as_i64(reader)?;
        let b4 = read_u8_as_i64(reader)?;
        let b5 = read_u8_as_i64(reader)?;
        let b6 = read_u8_as_i64(reader)?;
        (b0 & 0x03) << 48 | b1 << 40 | b2 << 32 | b3 << 24 | b4 << 16 | b5 << 8 | b6
    } else if b0 & 0x01 == 0 {
        let b1 = read_u8_as_i64(reader)?;
        let b2 = read_u8_as_i64(reader)?;
        let b3 = read_u8_as_i64(reader)?;
        let b4 = read_u8_as_i64(reader)?;
        let b5 = read_u8_as_i64(reader)?;
        let b6 = read_u8_as_i64(reader)?;
        let b7 = read_u8_as_i64(reader)?;
        b1 << 48 | b2 << 40 | b3 << 32 | b4 << 24 | b5 << 16 | b6 << 8 | b7
    } else {
        reader.read_i64::<BigEndian>()?
    };

    Ok(value)
}

pub fn write_ltf8<W>(writer: &mut W, value: Ltf8) -> io::Result<()>
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
    } else if value >> (40 - 5) == 0 {
        writer.write_u8((value >> 32 | 0xf0) as u8)?;
        writer.write_u8((value >> 24) as u8)?;
        writer.write_u8((value >> 16) as u8)?;
        writer.write_u8((value >> 8) as u8)?;
        writer.write_u8(value as u8)?;
    } else if value >> (48 - 6) == 0 {
        writer.write_u8((value >> 40 | 0xf8) as u8)?;
        writer.write_u8((value >> 32) as u8)?;
        writer.write_u8((value >> 24) as u8)?;
        writer.write_u8((value >> 16) as u8)?;
        writer.write_u8((value >> 8) as u8)?;
        writer.write_u8(value as u8)?;
    } else if value >> (56 - 7) == 0 {
        writer.write_u8((value >> 48 | 0xfc) as u8)?;
        writer.write_u8((value >> 40) as u8)?;
        writer.write_u8((value >> 32) as u8)?;
        writer.write_u8((value >> 24) as u8)?;
        writer.write_u8((value >> 16) as u8)?;
        writer.write_u8((value >> 8) as u8)?;
        writer.write_u8(value as u8)?;
    } else if value >> (64 - 8) == 0 {
        writer.write_u8((value >> 56 | 0xfe) as u8)?;
        writer.write_u8((value >> 48) as u8)?;
        writer.write_u8((value >> 40) as u8)?;
        writer.write_u8((value >> 32) as u8)?;
        writer.write_u8((value >> 24) as u8)?;
        writer.write_u8((value >> 16) as u8)?;
        writer.write_u8((value >> 8) as u8)?;
        writer.write_u8(value as u8)?;
    } else {
        writer.write_u8(0xff)?;
        writer.write_i64::<BigEndian>(value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_ltf8() -> io::Result<()> {
        let data = [0x00];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 0);

        let data = [0x55];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 85);

        let data = [0x80, 0xaa];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 170);

        let data = [0xc0, 0x55, 0xaa];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 21930);

        let data = [0xe0, 0x55, 0xaa, 0xcc];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 5614284);

        let data = [0xf0, 0x55, 0xaa, 0xcc, 0x33];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 1437256755);

        let data = [0xf8, 0x55, 0xaa, 0xcc, 0x33, 0xe3];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 367937729507);

        let data = [0xfc, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 94192058753820);

        let data = [0xfe, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 24113167040978160);

        let data = [0xff, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0, 0x0f];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, 6172970762490408975);

        let data = [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x56];
        let mut reader = &data[..];
        assert_eq!(read_ltf8(&mut reader)?, -170);

        Ok(())
    }

    #[test]
    fn test_write_ltf8() -> io::Result<()> {
        let mut buf = Vec::new();
        write_ltf8(&mut buf, 0)?;
        assert_eq!(buf, [0x00]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 85)?;
        assert_eq!(buf, [0x55]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 170)?;
        assert_eq!(buf, [0x80, 0xaa]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 21930)?;
        assert_eq!(buf, [0xc0, 0x55, 0xaa]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 5614284)?;
        assert_eq!(buf, [0xe0, 0x55, 0xaa, 0xcc]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 1437256755)?;
        assert_eq!(buf, [0xf0, 0x55, 0xaa, 0xcc, 0x33]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 367937729507)?;
        assert_eq!(buf, [0xf8, 0x55, 0xaa, 0xcc, 0x33, 0xe3]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 94192058753820)?;
        assert_eq!(buf, [0xfc, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 24113167040978160)?;
        assert_eq!(buf, [0xfe, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, 6172970762490408975)?;
        assert_eq!(buf, [0xff, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0, 0x0f]);

        let mut buf = Vec::new();
        write_ltf8(&mut buf, -170)?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x56]);

        Ok(())
    }
}
