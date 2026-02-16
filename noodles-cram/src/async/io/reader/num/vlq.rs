use std::num;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::io::reader::num::vlq::{zigzag_decode_i32, zigzag_decode_i64};

pub async fn read_uint7<R>(reader: &mut R) -> io::Result<u32>
where
    R: AsyncRead + Unpin,
{
    let mut n = 0u32;
    let mut count = 0u8;

    loop {
        let b = reader.read_u8().await.map(u32::from)?;

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

pub async fn read_uint7_as<R, N>(reader: &mut R) -> io::Result<N>
where
    R: AsyncRead + Unpin,
    N: TryFrom<u32, Error = num::TryFromIntError>,
{
    read_uint7(reader).await.and_then(|n| {
        n.try_into()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

/// Reads a signed integer using zigzag encoding (uint7 with zigzag decode).
pub async fn read_sint7<R>(reader: &mut R) -> io::Result<i32>
where
    R: AsyncRead + Unpin,
{
    let n = read_uint7(reader).await?;
    Ok(zigzag_decode_i32(n))
}

/// Reads a 64-bit unsigned VLQ integer.
pub async fn read_uint7_64<R>(reader: &mut R) -> io::Result<u64>
where
    R: AsyncRead + Unpin,
{
    let mut n: u64 = 0;
    let mut count = 0u8;

    loop {
        let b = reader.read_u8().await.map(u64::from)?;

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

#[allow(dead_code)]
pub async fn read_sint7_64<R>(reader: &mut R) -> io::Result<i64>
where
    R: AsyncRead + Unpin,
{
    let n = read_uint7_64(reader).await?;
    Ok(zigzag_decode_i64(n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_uint7() -> io::Result<()> {
        async fn t(mut data: &[u8], expected: u32) -> io::Result<()> {
            assert_eq!(read_uint7(&mut data).await?, expected);
            Ok(())
        }

        t(&[0x00], 0).await?;
        t(&[0x7f], 127).await?;
        t(&[0x81, 0x00], 128).await?;
        t(&[0xc0, 0x00], 8192).await?;
        t(&[0xff, 0x7f], 16383).await?;
        t(&[0x81, 0x80, 0x00], 16384).await?;

        Ok(())
    }

    #[tokio::test]
    async fn test_read_sint7() -> io::Result<()> {
        async fn t(mut data: &[u8], expected: i32) -> io::Result<()> {
            assert_eq!(read_sint7(&mut data).await?, expected);
            Ok(())
        }

        t(&[0x00], 0).await?;
        t(&[0x01], -1).await?;
        t(&[0x02], 1).await?;
        t(&[0x03], -2).await?;
        t(&[0x04], 2).await?;

        Ok(())
    }

    #[tokio::test]
    async fn test_read_uint7_64() -> io::Result<()> {
        async fn t(mut data: &[u8], expected: u64) -> io::Result<()> {
            assert_eq!(read_uint7_64(&mut data).await?, expected);
            Ok(())
        }

        t(&[0x00], 0).await?;
        t(&[0x7f], 127).await?;
        t(&[0x81, 0x00], 128).await?;

        Ok(())
    }

    #[tokio::test]
    async fn test_read_sint7_64() -> io::Result<()> {
        async fn t(mut data: &[u8], expected: i64) -> io::Result<()> {
            assert_eq!(read_sint7_64(&mut data).await?, expected);
            Ok(())
        }

        t(&[0x00], 0).await?;
        t(&[0x01], -1).await?;
        t(&[0x02], 1).await?;
        t(&[0x03], -2).await?;
        t(&[0x04], 2).await?;

        Ok(())
    }
}
