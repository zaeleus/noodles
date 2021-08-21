use std::io::{self, Write};

use byteorder::{BigEndian, WriteBytesExt};

use super::Ltf8;

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
