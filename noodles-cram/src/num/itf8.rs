use std::io::{self, Read};

use byteorder::ReadBytesExt;

fn read_u8_as_i32<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_u8().map(|b| b as i32)
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

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use super::*;

    #[test]
    fn test_read_itf8() -> io::Result<()> {
        let data = [0x00];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 0);

        let data = [0x87, 0x55];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 1877);

        let data = [0xc7, 0x55, 0x99];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 480665);

        let data = [0xe7, 0x55, 0x99, 0x66];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 123050342);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x02];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 1968805474);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x12];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 1968805474);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x22];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 1968805474);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x42];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 1968805474);

        let data = [0xf7, 0x55, 0x99, 0x66, 0x82];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, 1968805474);

        let data = [0xff, 0xff, 0xff, 0xff, 0x0f];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader)?;
        assert_eq!(i, -1);

        Ok(())
    }
}
