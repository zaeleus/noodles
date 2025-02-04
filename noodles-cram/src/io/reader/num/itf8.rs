use std::{
    io::{self, Read},
    num,
};

use byteorder::ReadBytesExt;
use bytes::Buf;

pub fn get_itf8<B>(src: &mut B) -> io::Result<i32>
where
    B: Buf,
{
    let b0 = src
        .try_get_u8()
        .map_err(|e| io::Error::new(io::ErrorKind::UnexpectedEof, e))?;

    // SAFETY: `leading_ones` is at max 4.
    let len = (b0 & 0xf0).leading_ones() as usize;

    if src.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let b0 = i32::from(b0);

    match len {
        0 => Ok(b0),
        1 => {
            let b1 = i32::from(src.get_u8());
            Ok((b0 & 0x7f) << 8 | b1)
        }
        2 => {
            let b1_2 = i32::from(src.get_u16());
            Ok((b0 & 0x3f) << 16 | b1_2)
        }
        3 => {
            let b1_2 = i32::from(src.get_u16());
            let b3 = i32::from(src.get_u8());
            Ok((b0 & 0x1f) << 24 | b1_2 << 8 | b3)
        }
        4 => {
            let b1_2 = i32::from(src.get_u16());
            let b3 = i32::from(src.get_u8());
            let b4 = i32::from(src.get_u8());
            Ok((b0 & 0x0f) << 28 | b1_2 << 12 | b3 << 4 | (b4 & 0x0f))
        }
        _ => unreachable!(),
    }
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
    fn test_get_itf8() -> io::Result<()> {
        fn t(mut data: &[u8], expected: i32) -> io::Result<()> {
            assert_eq!(get_itf8(&mut data)?, expected);
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
