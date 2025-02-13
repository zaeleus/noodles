use std::{
    io::{self, Read},
    num,
};

use byteorder::{BigEndian, ReadBytesExt};

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

fn read_u8_as_i64<R>(reader: &mut R) -> io::Result<i64>
where
    R: Read,
{
    reader.read_u8().map(i64::from)
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
