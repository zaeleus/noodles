use std::convert::TryFrom;

use noodles_bam as bam;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncReadExt};

use crate::{
    container::compression_header::{
        preservation_map::Key, PreservationMap, SubstitutionMatrix, TagIdsDictionary,
    },
    r#async::reader::num::read_itf8,
    record,
};

pub async fn read_preservation_map<R>(reader: &mut R) -> io::Result<PreservationMap>
where
    R: AsyncRead + Unpin,
{
    let data_len = read_itf8(reader).await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = vec![0; data_len];
    reader.read_exact(&mut buf).await?;

    let mut buf_reader = &buf[..];
    read_preservation_map_map(&mut buf_reader).await
}

async fn read_preservation_map_map<R>(reader: &mut R) -> io::Result<PreservationMap>
where
    R: AsyncRead + Unpin,
{
    let map_len = read_itf8(reader).await?;

    // § 8.4 "Compression header block" (2021-02-04): "The boolean values are optional, defaulting
    // to true when absent..."
    let mut read_names_included = true;
    let mut ap_data_series_delta = true;
    let mut reference_required = true;
    let mut substitution_matrix = None;
    let mut tag_ids_dictionary = None;

    for _ in 0..map_len {
        let key = read_key(reader).await?;

        match key {
            Key::ReadNamesIncluded => read_names_included = read_bool(reader).await?,
            Key::ApDataSeriesDelta => ap_data_series_delta = read_bool(reader).await?,
            Key::ReferenceRequired => reference_required = read_bool(reader).await?,
            Key::SubstitutionMatrix => {
                substitution_matrix = read_substitution_matrix(reader).await.map(Some)?;
            }
            Key::TagIdsDictionary => {
                tag_ids_dictionary = read_tag_ids_dictionary(reader).await.map(Some)?;
            }
        }
    }

    // § 8.4 "Compression header block" (2021-02-04): "SM and TD are mandatory."
    Ok(PreservationMap::new(
        read_names_included,
        ap_data_series_delta,
        reference_required,
        substitution_matrix.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing substitution matrix")
        })?,
        tag_ids_dictionary.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing tag IDs dictionary")
        })?,
    ))
}

async fn read_key<R>(reader: &mut R) -> io::Result<Key>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0; 2];
    reader.read_exact(&mut buf).await?;
    Key::try_from(&buf[..]).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

async fn read_bool<R>(reader: &mut R) -> io::Result<bool>
where
    R: AsyncRead + Unpin,
{
    match reader.read_u8().await {
        Ok(0) => Ok(false),
        Ok(1) => Ok(true),
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid bool value",
        )),
        Err(e) => Err(e),
    }
}

async fn read_substitution_matrix<R>(reader: &mut R) -> io::Result<SubstitutionMatrix>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0; 5];
    reader.read_exact(&mut buf[..]).await?;
    SubstitutionMatrix::try_from(&buf[..])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

async fn read_tag_ids_dictionary<R>(reader: &mut R) -> io::Result<TagIdsDictionary>
where
    R: AsyncRead + Unpin,
{
    let data_len = read_itf8(reader).await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = vec![0; data_len];
    reader.read_exact(&mut buf).await?;

    let mut buf_reader = &buf[..];
    let dictionary = read_raw_tag_ids_dictionary(&mut buf_reader).await?;

    Ok(TagIdsDictionary::from(dictionary))
}

async fn read_raw_tag_ids_dictionary<R>(reader: &mut R) -> io::Result<Vec<Vec<record::tag::Key>>>
where
    R: AsyncBufRead + Unpin,
{
    use bam::record::data::field::value::Type;

    use self::record::tag::Key;

    const NUL: u8 = 0x00;

    let mut dictionary = Vec::new();
    let mut keys_buf: Vec<u8> = Vec::new();

    loop {
        keys_buf.clear();

        match reader.read_until(NUL, &mut keys_buf).await {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(e),
        }

        let mut line = Vec::new();

        for chunk in keys_buf.chunks_exact(3) {
            let (t0, t1, ty) = (chunk[0], chunk[1], chunk[2]);

            let tag = [t0, t1];
            let ty =
                Type::try_from(ty).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let key = Key::new(tag, ty);

            line.push(key);
        }

        dictionary.push(line);
    }

    Ok(dictionary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_key() -> io::Result<()> {
        async fn t(mut reader: &[u8], expected: Key) -> io::Result<()> {
            let actual = read_key(&mut reader).await?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[b'R', b'N'], Key::ReadNamesIncluded).await?;
        t(&[b'A', b'P'], Key::ApDataSeriesDelta).await?;
        t(&[b'R', b'R'], Key::ReferenceRequired).await?;
        t(&[b'S', b'M'], Key::SubstitutionMatrix).await?;
        t(&[b'T', b'D'], Key::TagIdsDictionary).await?;

        let data = [b'Z', b'Z'];
        let mut reader = &data[..];
        assert!(matches!(
            read_key(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }

    #[tokio::test]
    async fn test_read_bool() -> io::Result<()> {
        let data = [0x00];
        let mut reader = &data[..];
        assert!(!read_bool(&mut reader).await?);

        let data = [0x01];
        let mut reader = &data[..];
        assert!(read_bool(&mut reader).await?);

        let data = [0x02];
        let mut reader = &data[..];
        assert!(matches!(
            read_bool(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[tokio::test]
    async fn test_read_subtitution_matrix() -> io::Result<()> {
        // [[C, G, T, N], [A, G, T, N], [A, C, T, N], [A, C, G, N], [A, C, G, T]]
        let data = [0x1b, 0x1b, 0x1b, 0x1b, 0x1b];
        let mut reader = &data[..];

        assert_eq!(
            read_substitution_matrix(&mut reader).await?,
            SubstitutionMatrix::default()
        );

        Ok(())
    }

    #[tokio::test]
    async fn test_read_tag_ids_dictionary() -> io::Result<()> {
        use noodles_bam::record::data::field::value::Type;

        use self::record::tag::Key;

        let data = [
            0x0b, // data_len = 11
            b'C', b'O', b'Z', // dictionary[0][0] = (b"CO", String)
            b'N', b'H', b'C', // dictionary[0][1] = (b"NH", UInt8)
            0x00, //
            b'P', b'G', b'Z', // dictionary[1][0] = (b"PG", String)
            0x00, //
        ];

        let mut reader = &data[..];
        let actual = read_tag_ids_dictionary(&mut reader).await?;

        let expected = TagIdsDictionary::from(vec![
            vec![
                Key::new([b'C', b'O'], Type::String),
                Key::new([b'N', b'H'], Type::UInt8),
            ],
            vec![Key::new([b'P', b'G'], Type::String)],
        ]);

        assert_eq!(actual, expected);

        Ok(())
    }
}
