use std::io::{self, Write};

pub fn write_uint7<W>(writer: &mut W, mut n: u32) -> io::Result<()>
where
    W: Write,
{
    let mut buf = [0; 5];

    let mut i = buf.len() - 1;
    let b = (n & 0x7f) as u8;
    buf[i] = b;
    n >>= 7;

    while n > 0 {
        i -= 1;
        let b = (n & 0x7f) as u8;
        buf[i] = b | 0x80;
        n >>= 7;
    }

    writer.write_all(&buf[i..])?;

    Ok(())
}

pub fn write_sint7<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    write_uint7(writer, zigzag_encode_i32(n))
}

pub fn write_uint7_64<W>(writer: &mut W, mut n: u64) -> io::Result<()>
where
    W: Write,
{
    let mut buf = [0; 10];

    let mut i = buf.len() - 1;
    let b = (n & 0x7f) as u8;
    buf[i] = b;
    n >>= 7;

    while n > 0 {
        i -= 1;
        let b = (n & 0x7f) as u8;
        buf[i] = b | 0x80;
        n >>= 7;
    }

    writer.write_all(&buf[i..])?;

    Ok(())
}

pub fn write_sint7_64<W>(writer: &mut W, n: i64) -> io::Result<()>
where
    W: Write,
{
    write_uint7_64(writer, zigzag_encode_i64(n))
}

fn zigzag_encode_i32(n: i32) -> u32 {
    ((n << 1) ^ (n >> 31)) as u32
}

fn zigzag_encode_i64(n: i64) -> u64 {
    ((n << 1) ^ (n >> 63)) as u64
}

/// Returns the encoded size of a uint7 value in bytes.
pub fn uint7_size_of(mut n: u32) -> usize {
    let mut size = 1;
    n >>= 7;
    while n > 0 {
        size += 1;
        n >>= 7;
    }
    size
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_uint7() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: u32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_uint7(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        // Examples from <https://en.wikipedia.org/wiki/Variable-length_quantity#Examples>.
        t(&mut buf, 0, &[0x00])?;
        t(&mut buf, 127, &[0x7f])?;
        t(&mut buf, 128, &[0x81, 0x00])?;
        t(&mut buf, 8192, &[0xc0, 0x00])?;
        t(&mut buf, 16383, &[0xff, 0x7f])?;
        t(&mut buf, 16384, &[0x81, 0x80, 0x00])?;
        t(&mut buf, 2097151, &[0xff, 0xff, 0x7f])?;
        t(&mut buf, 2097152, &[0x81, 0x80, 0x80, 0x00])?;
        t(&mut buf, 134217728, &[0xc0, 0x80, 0x80, 0x00])?;
        t(&mut buf, 268435455, &[0xff, 0xff, 0xff, 0x7f])?;

        Ok(())
    }

    #[test]
    fn test_write_sint7() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: i32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_sint7(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        // zigzag encoding: 0 -> 0, -1 -> 1, 1 -> 2, -2 -> 3, 2 -> 4, ...
        t(&mut buf, 0, &[0x00])?;
        t(&mut buf, -1, &[0x01])?;
        t(&mut buf, 1, &[0x02])?;
        t(&mut buf, -2, &[0x03])?;
        t(&mut buf, 2, &[0x04])?;

        Ok(())
    }

    #[test]
    fn test_write_uint7_64() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: u64, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_uint7_64(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &[0x00])?;
        t(&mut buf, 127, &[0x7f])?;
        t(&mut buf, 128, &[0x81, 0x00])?;

        Ok(())
    }

    #[test]
    fn test_write_sint7_64() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: i64, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_sint7_64(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &[0x00])?;
        t(&mut buf, -1, &[0x01])?;
        t(&mut buf, 1, &[0x02])?;
        t(&mut buf, -2, &[0x03])?;
        t(&mut buf, 2, &[0x04])?;

        Ok(())
    }

    #[test]
    fn test_zigzag_encode_i32() {
        assert_eq!(zigzag_encode_i32(0), 0);
        assert_eq!(zigzag_encode_i32(-1), 1);
        assert_eq!(zigzag_encode_i32(1), 2);
        assert_eq!(zigzag_encode_i32(-2), 3);
        assert_eq!(zigzag_encode_i32(2), 4);
        assert_eq!(zigzag_encode_i32(i32::MAX), 4294967294);
        assert_eq!(zigzag_encode_i32(i32::MIN), 4294967295);
    }

    #[test]
    fn test_zigzag_encode_i64() {
        assert_eq!(zigzag_encode_i64(0), 0);
        assert_eq!(zigzag_encode_i64(-1), 1);
        assert_eq!(zigzag_encode_i64(1), 2);
        assert_eq!(zigzag_encode_i64(-2), 3);
    }

    #[test]
    fn test_uint7_size_of() {
        assert_eq!(uint7_size_of(0), 1);
        assert_eq!(uint7_size_of(127), 1);
        assert_eq!(uint7_size_of(128), 2);
        assert_eq!(uint7_size_of(16383), 2);
        assert_eq!(uint7_size_of(16384), 3);
    }

    #[test]
    fn test_write_uint7_max() -> io::Result<()> {
        let mut buf = Vec::new();
        write_uint7(&mut buf, u32::MAX)?;
        assert_eq!(buf, &[0x8f, 0xff, 0xff, 0xff, 0x7f]);
        assert_eq!(uint7_size_of(u32::MAX), 5);
        Ok(())
    }

    #[test]
    fn test_write_uint7_64_large() -> io::Result<()> {
        let mut buf = Vec::new();
        write_uint7_64(&mut buf, u64::MAX)?;
        // u64::MAX requires 10 bytes in VLQ
        assert_eq!(buf.len(), 10);
        Ok(())
    }

    #[test]
    fn test_vlq_round_trip_u32() -> io::Result<()> {
        use crate::io::reader::num::read_uint7;

        let values = [0u32, 1, 127, 128, 255, 256, 16383, 16384, u32::MAX];
        for &val in &values {
            let mut buf = Vec::new();
            write_uint7(&mut buf, val)?;
            let mut src = &buf[..];
            let decoded = read_uint7(&mut src)?;
            assert_eq!(val, decoded, "round-trip failed for {val}");
        }
        Ok(())
    }

    #[test]
    fn test_vlq_round_trip_i32() -> io::Result<()> {
        use crate::io::reader::num::read_sint7;

        let values = [0i32, 1, -1, 2, -2, 127, -128, i32::MAX, i32::MIN];
        for &val in &values {
            let mut buf = Vec::new();
            write_sint7(&mut buf, val)?;
            let mut src = &buf[..];
            let decoded = read_sint7(&mut src)?;
            assert_eq!(val, decoded, "round-trip failed for {val}");
        }
        Ok(())
    }

    #[test]
    fn test_vlq_round_trip_u64() -> io::Result<()> {
        use crate::io::reader::num::read_uint7_64;

        let values = [0u64, 1, 127, 128, 16383, 16384, u32::MAX as u64, u64::MAX];
        for &val in &values {
            let mut buf = Vec::new();
            write_uint7_64(&mut buf, val)?;
            let mut src = &buf[..];
            let decoded = read_uint7_64(&mut src)?;
            assert_eq!(val, decoded, "round-trip failed for {val}");
        }
        Ok(())
    }

    #[test]
    fn test_vlq_round_trip_i64() -> io::Result<()> {
        use crate::io::reader::num::read_sint7_64;

        let values = [
            0i64,
            1,
            -1,
            127,
            -128,
            i32::MAX as i64,
            i32::MIN as i64,
            i64::MAX,
            i64::MIN,
        ];
        for &val in &values {
            let mut buf = Vec::new();
            write_sint7_64(&mut buf, val)?;
            let mut src = &buf[..];
            let decoded = read_sint7_64(&mut src)?;
            assert_eq!(val, decoded, "round-trip failed for {val}");
        }
        Ok(())
    }
}
