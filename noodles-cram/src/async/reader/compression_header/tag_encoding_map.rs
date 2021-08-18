use std::{collections::HashMap, convert::TryFrom};

use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{
    container::compression_header::{Encoding, TagEncodingMap},
    num::Itf8,
    r#async::reader::num::read_itf8,
};

use super::read_encoding;

pub async fn read_tag_encoding_map<R>(reader: &mut R) -> io::Result<TagEncodingMap>
where
    R: AsyncRead + Unpin,
{
    let data_len = read_itf8(reader).await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = vec![0; data_len];
    reader.read_exact(&mut buf).await?;

    let mut buf_reader = &buf[..];
    read_tag_encoding_map_map(&mut buf_reader)
        .await
        .map(TagEncodingMap::from)
}

async fn read_tag_encoding_map_map<R>(reader: &mut R) -> io::Result<HashMap<Itf8, Encoding>>
where
    R: AsyncRead + Unpin,
{
    let map_len = read_itf8(reader).await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut map = HashMap::with_capacity(map_len);

    for _ in 0..map_len {
        let key = read_itf8(reader).await?;
        let encoding = read_encoding(reader).await?;
        map.insert(key, encoding);
    }

    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_tag_encoding_map() -> io::Result<()> {
        let data = [
            0x0c, // data.len = 12
            0x01, // map.len = 1
            0xe0, 0x50, 0x47, 0x5a, // key[0] = PG:Z
            // value[0] = ByteArrayStop(b'\t', 5261146)
            0x05, 0x05, 0x09, 0xe0, 0x50, 0x47, 0x5a,
        ];

        let mut reader = &data[..];
        let actual = read_tag_encoding_map(&mut reader).await?;

        let expected = TagEncodingMap::from(
            vec![(5261146, Encoding::ByteArrayStop(b'\t', 5261146))]
                .into_iter()
                .collect::<HashMap<_, _>>(),
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
