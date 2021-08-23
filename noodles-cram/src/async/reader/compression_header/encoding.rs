use std::convert::TryFrom;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{data_container::compression_header::Encoding, r#async::reader::num::read_itf8};

pub async fn read_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: AsyncRead + Unpin,
{
    match read_itf8(reader).await? {
        0 => Ok(Encoding::Null),
        1 => read_external_encoding(reader).await,
        2 => unimplemented!("GOLOMB"),
        3 => read_huffman_encoding(reader).await,
        4 => todo!("BYTE_ARRAY_LEN"),
        5 => read_byte_array_stop_encoding(reader).await,
        6 => read_beta_encoding(reader).await,
        7 => read_subexp_encoding(reader).await,
        8 => unimplemented!("GOLOMB_RICE"),
        9 => read_gamma_encoding(reader).await,
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid encoding kind",
        )),
    }
}

async fn read_external_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: AsyncRead + Unpin,
{
    let args = read_args(reader).await?;
    let mut args_reader = &args[..];

    let block_content_id = read_itf8(&mut args_reader).await?;

    Ok(Encoding::External(block_content_id))
}

async fn read_huffman_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: AsyncRead + Unpin,
{
    let args = read_args(reader).await?;
    let mut args_reader = &args[..];

    let alphabet_len = read_itf8(&mut args_reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;
    let mut alphabet = Vec::with_capacity(alphabet_len);

    for _ in 0..alphabet_len {
        let symbol = read_itf8(&mut args_reader).await?;
        alphabet.push(symbol);
    }

    let bit_lens_len = read_itf8(&mut args_reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;
    let mut bit_lens = Vec::with_capacity(bit_lens_len);

    for _ in 0..bit_lens_len {
        let len = read_itf8(&mut args_reader).await?;
        bit_lens.push(len);
    }

    Ok(Encoding::Huffman(alphabet, bit_lens))
}

async fn read_byte_array_stop_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: AsyncRead + Unpin,
{
    let args = read_args(reader).await?;
    let mut args_reader = &args[..];

    let stop_byte = args_reader.read_u8().await?;
    let block_content_id = read_itf8(&mut args_reader).await?;

    Ok(Encoding::ByteArrayStop(stop_byte, block_content_id))
}

async fn read_beta_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: AsyncRead + Unpin,
{
    let args = read_args(reader).await?;
    let mut args_reader = &args[..];

    let offset = read_itf8(&mut args_reader).await?;
    let len = read_itf8(&mut args_reader).await?;

    Ok(Encoding::Beta(offset, len))
}

async fn read_subexp_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: AsyncRead + Unpin,
{
    let args = read_args(reader).await?;
    let mut args_reader = &args[..];

    let offset = read_itf8(&mut args_reader).await?;
    let k = read_itf8(&mut args_reader).await?;

    Ok(Encoding::Subexp(offset, k))
}

async fn read_gamma_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: AsyncRead + Unpin,
{
    let args = read_args(reader).await?;
    let mut args_reader = &args[..];

    let offset = read_itf8(&mut args_reader).await?;

    Ok(Encoding::Gamma(offset))
}

async fn read_args<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: AsyncRead + Unpin,
{
    let len = read_itf8(reader).await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;
    let mut buf = vec![0; len];
    reader.read_exact(&mut buf).await?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_null_encoding() -> io::Result<()> {
        let data = [
            0x00, // codec ID = NULL
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader).await?;
        assert_eq!(encoding, Encoding::Null);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_external_encoding() -> io::Result<()> {
        let data = [
            0x01, // codec ID = EXTERNAL
            0x01, // args.len = 1
            0x05, // block content ID = 5
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader).await?;
        assert_eq!(encoding, Encoding::External(5));

        Ok(())
    }

    #[tokio::test]
    async fn test_read_huffman_encoding() -> io::Result<()> {
        let data = [
            0x03, // codec ID = HUFFMAN
            0x04, // args.len = 4
            0x01, // alphabet.len = 1
            0x41, // alphabet[0] = b'A'
            0x01, // bit_lens.len = 1
            0x00, // bit_lens[0] = 0
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader).await?;
        assert_eq!(encoding, Encoding::Huffman(vec![0x41], vec![0]));

        Ok(())
    }

    /* #[tokio::test]
    async fn test_read_byte_array_len_encoding() -> io::Result<()> {
        let data = [
            0x04, // codec ID = BYTE_ARRAY_LEN
            0x06, // args.len = 6
            0x01, // codec ID = EXTERNAL
            0x01, // args.len = 1
            0x0d, // block content ID = 13
            0x01, // codec ID = EXTERNAL
            0x01, // args.len = 1
            0x15, // block content ID = 21
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader).await?;

        assert_eq!(
            encoding,
            Encoding::ByteArrayLen(
                Box::new(Encoding::External(13)),
                Box::new(Encoding::External(21))
            )
        );

        Ok(())
    } */

    #[tokio::test]
    async fn test_read_byte_array_stop_encoding() -> io::Result<()> {
        let data = [
            0x05, // codec ID = BYTE_ARRAY_STOP
            0x02, // args.len = 2
            0x00, // stop = NUL
            0x08, // block content ID = 8
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader).await?;
        assert_eq!(encoding, Encoding::ByteArrayStop(0x00, 8));

        Ok(())
    }

    #[tokio::test]
    async fn test_read_beta_encoding() -> io::Result<()> {
        let data = [
            0x06, // codec ID = BETA
            0x02, // args.len = 2
            0x00, // offset = 0
            0x08, // len = 8
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader).await?;
        assert_eq!(encoding, Encoding::Beta(0, 8));

        Ok(())
    }

    #[tokio::test]
    async fn test_read_subexp_encoding() -> io::Result<()> {
        let data = [
            0x07, // codec ID = SUBEXP
            0x02, // args.len = 2
            0x00, // offset = 0
            0x01, // k = 1
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader).await?;
        assert_eq!(encoding, Encoding::Subexp(0, 1));

        Ok(())
    }

    #[tokio::test]
    async fn test_read_gamma_encoding() -> io::Result<()> {
        let data = [
            0x09, // codec ID = GAMMA
            0x01, // args.len
            0x01, // offset
        ];
        let mut reader = &data[..];

        let encoding = read_encoding(&mut reader).await?;
        assert_eq!(encoding, Encoding::Gamma(1));

        Ok(())
    }
}
