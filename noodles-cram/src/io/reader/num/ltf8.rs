use std::{
    io::{self, Read},
    num,
};

use super::{read_u8, read_u16_be, read_u24_be, read_u32_be};

pub fn read_ltf8<R>(reader: &mut R) -> io::Result<i64>
where
    R: Read,
{
    let b0 = read_u8(reader).map(u64::from)?;

    let n = if b0 & 0x80 == 0 {
        b0
    } else if b0 & 0x40 == 0 {
        let b1 = read_u8(reader).map(u64::from)?;
        (b0 & 0x7f) << 8 | b1
    } else if b0 & 0x20 == 0 {
        let b1_2 = read_u16_be(reader).map(u64::from)?;
        (b0 & 0x3f) << 16 | b1_2
    } else if b0 & 0x10 == 0 {
        let b1_3 = read_u24_be(reader).map(u64::from)?;
        (b0 & 0x1f) << 24 | b1_3
    } else if b0 & 0x08 == 0 {
        let b1_4 = read_u32_be(reader).map(u64::from)?;
        (b0 & 0x0f) << 32 | b1_4
    } else if b0 & 0x04 == 0 {
        let b1_5 = read_u40_be(reader)?;
        (b0 & 0x07) << 40 | b1_5
    } else if b0 & 0x02 == 0 {
        let b1_6 = read_u48_be(reader)?;
        (b0 & 0x03) << 48 | b1_6
    } else if b0 & 0x01 == 0 {
        read_u56_be(reader)?
    } else {
        read_u64_be(reader)?
    };

    Ok(n as i64)
}

pub fn read_ltf8_as<R, N>(reader: &mut R) -> io::Result<N>
where
    R: Read,
    N: TryFrom<i64, Error = num::TryFromIntError>,
{
    read_ltf8(reader).and_then(|n| {
        n.try_into()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_u40_be<R>(reader: &mut R) -> io::Result<u64>
where
    R: Read,
{
    let mut buf = [0; 8];
    reader.read_exact(&mut buf[3..])?;
    Ok(u64::from_be_bytes(buf))
}

fn read_u48_be<R>(reader: &mut R) -> io::Result<u64>
where
    R: Read,
{
    let mut buf = [0; 8];
    reader.read_exact(&mut buf[2..])?;
    Ok(u64::from_be_bytes(buf))
}

fn read_u56_be<R>(reader: &mut R) -> io::Result<u64>
where
    R: Read,
{
    let mut buf = [0; 8];
    reader.read_exact(&mut buf[1..])?;
    Ok(u64::from_be_bytes(buf))
}

fn read_u64_be<R>(reader: &mut R) -> io::Result<u64>
where
    R: Read,
{
    let mut buf = [0; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_be_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_ltf8() -> io::Result<()> {
        fn t(mut data: &[u8], expected: i64) -> io::Result<()> {
            assert_eq!(read_ltf8(&mut data)?, expected);
            Ok(())
        }

        t(&[0x00], 0)?;
        t(&[0x55], 85)?;
        t(&[0x80, 0xaa], 170)?;
        t(&[0xc0, 0x55, 0xaa], 21930)?;
        t(&[0xe0, 0x55, 0xaa, 0xcc], 5614284)?;
        t(&[0xf0, 0x55, 0xaa, 0xcc, 0x33], 1437256755)?;
        t(&[0xf8, 0x55, 0xaa, 0xcc, 0x33, 0xe3], 367937729507)?;
        t(&[0xfc, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c], 94192058753820)?;
        t(
            &[0xfe, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0],
            24113167040978160,
        )?;
        t(
            &[0xff, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0, 0x0f],
            6172970762490408975,
        )?;
        t(
            &[0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x56],
            -170,
        )?;

        Ok(())
    }
}
