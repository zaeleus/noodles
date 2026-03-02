use std::{
    io::{self, Read},
    num,
};

use super::read_u8;

pub fn read_uint7<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut n = 0u32;
    let mut count = 0u8;

    loop {
        let b = read_u8(reader).map(u32::from)?;

        count += 1;
        if count > 5 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "VLQ integer overflow",
            ));
        }

        n <<= 7;
        n |= b & 0x7f;

        if b & 0x80 == 0 {
            break;
        }
    }

    Ok(n)
}

pub fn read_uint7_as<R, N>(reader: &mut R) -> io::Result<N>
where
    R: Read,
    N: TryFrom<u32, Error = num::TryFromIntError>,
{
    read_uint7(reader).and_then(|n| {
        n.try_into()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

/// Reads a signed integer using zigzag encoding (uint7 with zigzag decode).
///
/// Used in CRAM 4.0 for signed integer fields.
pub fn read_sint7<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    let n = read_uint7(reader)?;
    Ok(zigzag_decode_i32(n))
}

/// Reads a 64-bit unsigned VLQ integer.
///
/// Used in CRAM 4.0 for long fields (e.g., record counter, base count).
pub fn read_uint7_64<R>(reader: &mut R) -> io::Result<u64>
where
    R: Read,
{
    let mut n: u64 = 0;
    let mut count = 0u8;

    loop {
        let b = read_u8(reader).map(u64::from)?;

        count += 1;
        if count > 10 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "VLQ integer overflow",
            ));
        }

        n <<= 7;
        n |= b & 0x7f;

        if b & 0x80 == 0 {
            break;
        }
    }

    Ok(n)
}

/// Reads a signed 64-bit integer using zigzag encoding (uint7_64 with zigzag decode).
///
/// Used in CRAM 4.0 for signed 64-bit fields.
pub fn read_sint7_64<R>(reader: &mut R) -> io::Result<i64>
where
    R: Read,
{
    let n = read_uint7_64(reader)?;
    Ok(zigzag_decode_i64(n))
}

pub(crate) fn zigzag_decode_i32(n: u32) -> i32 {
    ((n >> 1) as i32) ^ -((n & 1) as i32)
}

pub(crate) fn zigzag_decode_i64(n: u64) -> i64 {
    ((n >> 1) as i64) ^ -((n & 1) as i64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_uint7() -> io::Result<()> {
        fn t(mut data: &[u8], expected: u32) -> io::Result<()> {
            assert_eq!(read_uint7(&mut data)?, expected);
            Ok(())
        }

        // Examples from <https://en.wikipedia.org/wiki/Variable-length_quantity#Examples>.
        t(&[0x00], 0)?;
        t(&[0x7f], 127)?;
        t(&[0x81, 0x00], 128)?;
        t(&[0xc0, 0x00], 8192)?;
        t(&[0xff, 0x7f], 16383)?;
        t(&[0x81, 0x80, 0x00], 16384)?;
        t(&[0xff, 0xff, 0x7f], 2097151)?;
        t(&[0x81, 0x80, 0x80, 0x00], 2097152)?;
        t(&[0xc0, 0x80, 0x80, 0x00], 134217728)?;
        t(&[0xff, 0xff, 0xff, 0x7f], 268435455)?;

        Ok(())
    }

    #[test]
    fn test_read_sint7() -> io::Result<()> {
        fn t(mut data: &[u8], expected: i32) -> io::Result<()> {
            assert_eq!(read_sint7(&mut data)?, expected);
            Ok(())
        }

        // zigzag encoding: 0 -> 0, 1 -> -1, 2 -> 1, 3 -> -2, 4 -> 2, ...
        t(&[0x00], 0)?;
        t(&[0x01], -1)?;
        t(&[0x02], 1)?;
        t(&[0x03], -2)?;
        t(&[0x04], 2)?;

        Ok(())
    }

    #[test]
    fn test_read_uint7_64() -> io::Result<()> {
        fn t(mut data: &[u8], expected: u64) -> io::Result<()> {
            assert_eq!(read_uint7_64(&mut data)?, expected);
            Ok(())
        }

        t(&[0x00], 0)?;
        t(&[0x7f], 127)?;
        t(&[0x81, 0x00], 128)?;

        Ok(())
    }

    #[test]
    fn test_read_sint7_64() -> io::Result<()> {
        fn t(mut data: &[u8], expected: i64) -> io::Result<()> {
            assert_eq!(read_sint7_64(&mut data)?, expected);
            Ok(())
        }

        t(&[0x00], 0)?;
        t(&[0x01], -1)?;
        t(&[0x02], 1)?;
        t(&[0x03], -2)?;
        t(&[0x04], 2)?;

        Ok(())
    }

    #[test]
    fn test_read_uint7_overflow() {
        // 6 continuation bytes (all with high bit set) should overflow a u32 VLQ
        let data: &[u8] = &[0x80, 0x80, 0x80, 0x80, 0x80, 0x00];
        let result = read_uint7(&mut &data[..]);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn test_read_uint7_64_overflow() {
        // 11 continuation bytes should overflow a u64 VLQ
        let data: &[u8] = &[
            0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00,
        ];
        let result = read_uint7_64(&mut &data[..]);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn test_zigzag_decode_i32() {
        assert_eq!(zigzag_decode_i32(0), 0);
        assert_eq!(zigzag_decode_i32(1), -1);
        assert_eq!(zigzag_decode_i32(2), 1);
        assert_eq!(zigzag_decode_i32(3), -2);
        assert_eq!(zigzag_decode_i32(4), 2);
        assert_eq!(zigzag_decode_i32(4294967294), 2147483647); // i32::MAX
        assert_eq!(zigzag_decode_i32(4294967295), i32::MIN); // i32::MIN
    }

    #[test]
    fn test_zigzag_decode_i64() {
        assert_eq!(zigzag_decode_i64(0), 0);
        assert_eq!(zigzag_decode_i64(1), -1);
        assert_eq!(zigzag_decode_i64(2), 1);
        assert_eq!(zigzag_decode_i64(3), -2);
    }
}
