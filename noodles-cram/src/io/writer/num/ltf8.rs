use std::io::{self, Write};

use byteorder::{BigEndian, WriteBytesExt};

pub fn write_ltf8<W>(writer: &mut W, n: i64) -> io::Result<()>
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
    } else if n >> (40 - 5) == 0 {
        let buf = [
            ((n >> 32) | 0xf0) as u8,
            (n >> 24) as u8,
            (n >> 16) as u8,
            (n >> 8) as u8,
            n as u8,
        ];

        writer.write_all(&buf)?;
    } else if n >> (48 - 6) == 0 {
        let buf = [
            ((n >> 40) | 0xf8) as u8,
            (n >> 32) as u8,
            (n >> 24) as u8,
            (n >> 16) as u8,
            (n >> 8) as u8,
            n as u8,
        ];

        writer.write_all(&buf)?;
    } else if n >> (56 - 7) == 0 {
        let buf = [
            ((n >> 48) | 0xfc) as u8,
            (n >> 40) as u8,
            (n >> 32) as u8,
            (n >> 24) as u8,
            (n >> 16) as u8,
            (n >> 8) as u8,
            n as u8,
        ];

        writer.write_all(&buf)?;
    } else if n >> (64 - 8) == 0 {
        let buf = [
            ((n >> 56) | 0xfe) as u8,
            (n >> 48) as u8,
            (n >> 40) as u8,
            (n >> 32) as u8,
            (n >> 24) as u8,
            (n >> 16) as u8,
            (n >> 8) as u8,
            n as u8,
        ];

        writer.write_all(&buf)?;
    } else {
        writer.write_u8(0xff)?;
        writer.write_i64::<BigEndian>(n)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_ltf8() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, value: i64, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_ltf8(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &[0x00])?;
        t(&mut buf, 85, &[0x55])?;
        t(&mut buf, 170, &[0x80, 0xaa])?;
        t(&mut buf, 21930, &[0xc0, 0x55, 0xaa])?;
        t(&mut buf, 5614284, &[0xe0, 0x55, 0xaa, 0xcc])?;
        t(&mut buf, 1437256755, &[0xf0, 0x55, 0xaa, 0xcc, 0x33])?;
        t(
            &mut buf,
            367937729507,
            &[0xf8, 0x55, 0xaa, 0xcc, 0x33, 0xe3],
        )?;
        t(
            &mut buf,
            94192058753820,
            &[0xfc, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c],
        )?;
        t(
            &mut buf,
            24113167040978160,
            &[0xfe, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0],
        )?;
        t(
            &mut buf,
            6172970762490408975,
            &[0xff, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0, 0x0f],
        )?;
        t(
            &mut buf,
            -170,
            &[0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x56],
        )?;

        Ok(())
    }
}
