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

    let n = (b0 & 0xf0).count_ones();

    let value = match n {
        0 => b0,
        1 => {
            let b1 = read_u8_as_i32(reader)?;
            (b0 & 0x7f) << 8 | b1
        }
        2 => {
            let b1 = read_u8_as_i32(reader)?;
            let b2 = read_u8_as_i32(reader)?;
            (b0 & 0x3f) << 16 | b1 << 8 | b2
        }
        3 => {
            let b1 = read_u8_as_i32(reader)?;
            let b2 = read_u8_as_i32(reader)?;
            let b3 = read_u8_as_i32(reader)?;
            (b0 & 0x1f) << 24 | b1 << 16 | b2 << 8 | b3
        }
        4 => {
            let b1 = read_u8_as_i32(reader)?;
            let b2 = read_u8_as_i32(reader)?;
            let b3 = read_u8_as_i32(reader)?;
            let b4 = read_u8_as_i32(reader)?;
            (b0 & 0x07) << 28 | b1 << 20 | b2 << 12 | b3 << 4 | b4 & 0x07
        }
        _ => unreachable!(),
    };

    Ok(value)
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;

    use super::*;

    #[test]
    fn test_read_itf8() {
        let data = [0x00];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader).unwrap();
        assert_eq!(i, 0);

        let data = [0x87, 0x55];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader).unwrap();
        assert_eq!(i, 1877);

        let data = [0xc7, 0x55, 0x99];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader).unwrap();
        assert_eq!(i, 480665);

        let data = [0xe7, 0x55, 0x99, 0x66];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader).unwrap();
        assert_eq!(i, 123050342);

        let data = [0xf7, 0x55, 0x99, 0x66, 0xd2];
        let mut reader = BufReader::new(&data[..]);
        let i = read_itf8(&mut reader).unwrap();
        assert_eq!(i, 1968805474);
    }
}
