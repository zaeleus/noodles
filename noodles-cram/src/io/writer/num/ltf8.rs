use std::io::{self, Write};

pub fn write_ltf8<W>(writer: &mut W, n: i64) -> io::Result<()>
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
    } else if n >> (40 - 5) == 0 {
        let m = (n as u64) | (0xf0 << 32);
        let buf = m.to_be_bytes();
        writer.write_all(&buf[3..])?;
    } else if n >> (48 - 6) == 0 {
        let m = (n as u64) | (0xf8 << 40);
        let buf = m.to_be_bytes();
        writer.write_all(&buf[2..])?;
    } else if n >> (56 - 7) == 0 {
        let m = (n as u64) | (0xfc << 48);
        let buf = m.to_be_bytes();
        writer.write_all(&buf[1..])?;
    } else if n >> (64 - 8) == 0 {
        let m = (n as u64) | (0xfe << 56);
        let buf = m.to_be_bytes();
        writer.write_all(&buf)?;
    } else {
        writer.write_all(&[0xff])?;
        write_i64_be(writer, n)?;
    }

    Ok(())
}

fn write_i64_be<W>(writer: &mut W, n: i64) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_be_bytes();
    writer.write_all(&buf)
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
