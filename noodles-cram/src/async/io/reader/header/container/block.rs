use async_compression::tokio::bufread::GzipDecoder;
use tokio::io::{self, AsyncRead, AsyncReadExt, BufReader};

use crate::{
    container::block::{CompressionMethod, ContentType},
    r#async::io::reader::num::{read_itf8, read_itf8_as},
};

pub(super) async fn read_block<R>(reader: &mut R) -> io::Result<Box<dyn AsyncRead + Unpin + '_>>
where
    R: AsyncRead + Unpin,
{
    let compression_method = read_compression_method(reader).await?;

    let content_type = read_content_type(reader).await?;
    validate_content_type(content_type)?;

    let _content_type_id = read_itf8(reader).await?;
    let compressed_size = read_itf8_as(reader).await?;
    let uncompressed_size = read_itf8_as(reader).await?;

    let reader: Box<dyn AsyncRead + Unpin + '_> = match compression_method {
        CompressionMethod::None => Box::new(reader.take(uncompressed_size)),
        CompressionMethod::Gzip => Box::new(GzipDecoder::new(BufReader::new(
            reader.take(compressed_size),
        ))),
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid block compression method: expected {:?} or {:?}, got {:?}",
                    CompressionMethod::None,
                    CompressionMethod::Gzip,
                    compression_method,
                ),
            ))
        }
    };

    // Skip CRC32.

    Ok(reader)
}

async fn read_compression_method<R>(reader: &mut R) -> io::Result<CompressionMethod>
where
    R: AsyncRead + Unpin,
{
    match reader.read_u8().await? {
        0 => Ok(CompressionMethod::None),
        1 => Ok(CompressionMethod::Gzip),
        2 => Ok(CompressionMethod::Bzip2),
        3 => Ok(CompressionMethod::Lzma),
        4 => Ok(CompressionMethod::Rans4x8),
        5 => Ok(CompressionMethod::RansNx16),
        6 => Ok(CompressionMethod::AdaptiveArithmeticCoding),
        7 => Ok(CompressionMethod::Fqzcomp),
        8 => Ok(CompressionMethod::NameTokenizer),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid compression method",
        )),
    }
}

async fn read_content_type<R>(reader: &mut R) -> io::Result<ContentType>
where
    R: AsyncRead + Unpin,
{
    match reader.read_u8().await? {
        0 => Ok(ContentType::FileHeader),
        1 => Ok(ContentType::CompressionHeader),
        2 => Ok(ContentType::SliceHeader),
        3 => Ok(ContentType::Reserved),
        4 => Ok(ContentType::ExternalData),
        5 => Ok(ContentType::CoreData),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid content type",
        )),
    }
}

fn validate_content_type(actual: ContentType) -> io::Result<()> {
    const EXPECTED: ContentType = ContentType::FileHeader;

    if actual == EXPECTED {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid block content type: expected {EXPECTED:?}, got {actual:?}"),
        ))
    }
}