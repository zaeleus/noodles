use std::{
    io::{self, Read},
    num,
};

pub fn read_itf8<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    let b0 = read_u8(reader).map(u32::from)?;

    let n = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        let b1 = read_u8(reader).map(u32::from)?;
        (b0 & 0x7f) << 8 | b1
    } else if b0 & 0x20 == 0 {
        let b1_2 = read_u16_be(reader).map(u32::from)?;
        (b0 & 0x3f) << 16 | b1_2
    } else if b0 & 0x10 == 0 {
        let b1_3 = read_u24_be(reader)?;
        (b0 & 0x1f) << 24 | b1_3
    } else {
        let b1_4 = read_u32_be(reader)?;
        (b0 & 0x0f) << 28 | (b1_4 & 0xffffff0f) >> 4 | b1_4 & 0x0f
    };

    Ok(n as i32)
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

fn read_u8<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u16_be<R>(reader: &mut R) -> io::Result<u16>
where
    R: Read,
{
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_be_bytes(buf))
}

fn read_u24_be<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; 3];
    reader.read_exact(&mut buf)?;
    Ok(u32::from(buf[0]) << 16 | u32::from(buf[1]) << 8 | u32::from(buf[2]))
}

fn read_u32_be<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
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
