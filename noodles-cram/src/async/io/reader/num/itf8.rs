use tokio::io::{self, AsyncRead, AsyncReadExt};

pub async fn read_itf8<R>(reader: &mut R) -> io::Result<i32>
where
    R: AsyncRead + Unpin,
{
    let b0 = read_u8_as_i32(reader).await?;

    if b0 & 0x80 == 0 {
        Ok(b0)
    } else if b0 & 0x40 == 0 {
        let b1 = read_u8_as_i32(reader).await?;
        Ok(((b0 & 0x7f) << 8) | b1)
    } else if b0 & 0x20 == 0 {
        let b1 = read_u8_as_i32(reader).await?;
        let b2 = read_u8_as_i32(reader).await?;
        Ok(((b0 & 0x3f) << 16) | (b1 << 8) | b2)
    } else if b0 & 0x10 == 0 {
        let b1 = read_u8_as_i32(reader).await?;
        let b2 = read_u8_as_i32(reader).await?;
        let b3 = read_u8_as_i32(reader).await?;
        Ok(((b0 & 0x1f) << 24) | (b1 << 16) | (b2 << 8) | b3)
    } else {
        let b1 = read_u8_as_i32(reader).await?;
        let b2 = read_u8_as_i32(reader).await?;
        let b3 = read_u8_as_i32(reader).await?;
        let b4 = read_u8_as_i32(reader).await?;
        Ok(((b0 & 0x0f) << 28) | (b1 << 20) | (b2 << 12) | (b3 << 4) | b4 & 0x0f)
    }
}

async fn read_u8_as_i32<R>(reader: &mut R) -> io::Result<i32>
where
    R: AsyncRead + Unpin,
{
    reader.read_u8().await.map(i32::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_itf8() -> io::Result<()> {
        async fn t(mut reader: &[u8], expected: i32) -> io::Result<()> {
            let actual = read_itf8(&mut reader).await?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[0x00], 0).await?;
        t(&[0x7f], 127).await?;

        t(&[0x87, 0x55], 1877).await?;
        t(&[0xc7, 0x55, 0x99], 480665).await?;
        t(&[0xe7, 0x55, 0x99, 0x66], 123050342).await?;

        t(&[0xf7, 0x55, 0x99, 0x66, 0x02], 1968805474).await?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x12], 1968805474).await?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x22], 1968805474).await?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x42], 1968805474).await?;
        t(&[0xf7, 0x55, 0x99, 0x66, 0x82], 1968805474).await?;

        t(&[0xff, 0xff, 0xff, 0xff, 0x0f], -1).await?;

        Ok(())
    }
}
