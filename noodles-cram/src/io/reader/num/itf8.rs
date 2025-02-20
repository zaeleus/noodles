use std::{
    io::{self, Read},
    num,
};

use byteorder::ReadBytesExt;

pub fn read_itf8<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    let b0 = read_u8_as_i32(reader)?;

    let value = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        let b1 = read_u8_as_i32(reader)?;
        ((b0 & 0x7f) << 8) | b1
    } else if b0 & 0x20 == 0 {
        let b1 = read_u8_as_i32(reader)?;
        let b2 = read_u8_as_i32(reader)?;
        ((b0 & 0x3f) << 16) | (b1 << 8) | b2
    } else if b0 & 0x10 == 0 {
        let b1 = read_u8_as_i32(reader)?;
        let b2 = read_u8_as_i32(reader)?;
        let b3 = read_u8_as_i32(reader)?;
        ((b0 & 0x1f) << 24) | (b1 << 16) | (b2 << 8) | b3
    } else {
        let b1 = read_u8_as_i32(reader)?;
        let b2 = read_u8_as_i32(reader)?;
        let b3 = read_u8_as_i32(reader)?;
        let b4 = read_u8_as_i32(reader)?;
        ((b0 & 0x0f) << 28) | (b1 << 20) | (b2 << 12) | (b3 << 4) | b4 & 0x0f
    };

    Ok(value)
}

pub fn read_itf8_as<R, N>(reader: &mut R) -> io::Result<N>
where
    R: Read,
    N: TryFrom<i32, Error = num::TryFromIntError>,
{
    read_itf8(reader).and_then(|n| {
        n.try_into()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_u8_as_i32<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_u8().map(i32::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_itf8() -> io::Result<()> {
        fn t(mut data: &[u8], expected: i32) -> io::Result<()> {
            assert_eq!(read_itf8(&mut data)?, expected);
            Ok(())
        }

        t(&[0x00], 0)?;
        t(&[0x87, 0x55], 1877)?;
        t(&[0xc7, 0x55, 0x99], 480665)?;
        t(&[0xe7, 0x55, 0x99, 0x66], 123050342)?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x02], 1968805474)?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x12], 1968805474)?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x22], 1968805474)?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x42], 1968805474)?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x82], 1968805474)?;
        t(&[0xff, 0xff, 0xff, 0xff, 0x0f], -1)?;

        Ok(())
    }
}
