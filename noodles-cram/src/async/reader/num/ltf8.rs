use tokio::io::{self, AsyncRead, AsyncReadExt};

pub async fn read_ltf8<R>(reader: &mut R) -> io::Result<i64>
where
    R: AsyncRead + Unpin,
{
    let b0 = read_u8_as_i64(reader).await?;

    if b0 & 0x80 == 0 {
        Ok(b0)
    } else if b0 & 0x40 == 0 {
        let b1 = read_u8_as_i64(reader).await?;
        Ok((b0 & 0x7f) << 8 | b1)
    } else if b0 & 0x20 == 0 {
        let b1 = read_u8_as_i64(reader).await?;
        let b2 = read_u8_as_i64(reader).await?;
        Ok((b0 & 0x3f) << 16 | b1 << 8 | b2)
    } else if b0 & 0x10 == 0 {
        let b1 = read_u8_as_i64(reader).await?;
        let b2 = read_u8_as_i64(reader).await?;
        let b3 = read_u8_as_i64(reader).await?;
        Ok((b0 & 0x1f) << 24 | b1 << 16 | b2 << 8 | b3)
    } else if b0 & 0x08 == 0 {
        let b1 = read_u8_as_i64(reader).await?;
        let b2 = read_u8_as_i64(reader).await?;
        let b3 = read_u8_as_i64(reader).await?;
        let b4 = read_u8_as_i64(reader).await?;
        Ok((b0 & 0x0f) << 32 | b1 << 24 | b2 << 16 | b3 << 8 | b4)
    } else if b0 & 0x04 == 0 {
        let b1 = read_u8_as_i64(reader).await?;
        let b2 = read_u8_as_i64(reader).await?;
        let b3 = read_u8_as_i64(reader).await?;
        let b4 = read_u8_as_i64(reader).await?;
        let b5 = read_u8_as_i64(reader).await?;
        Ok((b0 & 0x07) << 40 | b1 << 32 | b2 << 24 | b3 << 16 | b4 << 8 | b5)
    } else if b0 & 0x02 == 0 {
        let b1 = read_u8_as_i64(reader).await?;
        let b2 = read_u8_as_i64(reader).await?;
        let b3 = read_u8_as_i64(reader).await?;
        let b4 = read_u8_as_i64(reader).await?;
        let b5 = read_u8_as_i64(reader).await?;
        let b6 = read_u8_as_i64(reader).await?;
        Ok((b0 & 0x03) << 48 | b1 << 40 | b2 << 32 | b3 << 24 | b4 << 16 | b5 << 8 | b6)
    } else if b0 & 0x01 == 0 {
        let b1 = read_u8_as_i64(reader).await?;
        let b2 = read_u8_as_i64(reader).await?;
        let b3 = read_u8_as_i64(reader).await?;
        let b4 = read_u8_as_i64(reader).await?;
        let b5 = read_u8_as_i64(reader).await?;
        let b6 = read_u8_as_i64(reader).await?;
        let b7 = read_u8_as_i64(reader).await?;
        Ok(b1 << 48 | b2 << 40 | b3 << 32 | b4 << 24 | b5 << 16 | b6 << 8 | b7)
    } else {
        reader.read_i64().await
    }
}

async fn read_u8_as_i64<R>(reader: &mut R) -> io::Result<i64>
where
    R: AsyncRead + Unpin,
{
    reader.read_u8().await.map(i64::from)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_ltf8() -> io::Result<()> {
        async fn t(mut reader: &[u8], expected: i64) -> io::Result<()> {
            let actual = read_ltf8(&mut reader).await?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[0x00], 0).await?;
        t(&[0x55], 85).await?;
        t(&[0x80, 0xaa], 170).await?;
        t(&[0xc0, 0x55, 0xaa], 21930).await?;
        t(&[0xe0, 0x55, 0xaa, 0xcc], 5614284).await?;
        t(&[0xf0, 0x55, 0xaa, 0xcc, 0x33], 1437256755).await?;
        t(&[0xf8, 0x55, 0xaa, 0xcc, 0x33, 0xe3], 367937729507).await?;
        t(&[0xfc, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c], 94192058753820).await?;

        t(
            &[0xfe, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0],
            24113167040978160,
        )
        .await?;

        t(
            &[0xff, 0x55, 0xaa, 0xcc, 0x33, 0xe3, 0x1c, 0xf0, 0x0f],
            6172970762490408975,
        )
        .await?;

        t(
            &[0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x56],
            -170,
        )
        .await?;

        Ok(())
    }
}
