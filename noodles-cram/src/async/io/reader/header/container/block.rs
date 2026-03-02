use std::{
    pin::Pin,
    task::{Context, Poll},
};

use async_compression::tokio::bufread::GzipDecoder;
use pin_project_lite::pin_project;
use tokio::io::{self, AsyncRead, AsyncReadExt, BufReader, ReadBuf, Take};

use crate::{
    r#async::io::reader::num::{read_signed_int, read_unsigned_int_as},
    container::block::{CompressionMethod, ContentType},
    file_definition::Version,
};

pub(super) async fn read_block<R>(reader: &mut R, version: Version) -> io::Result<Decoder<&mut R>>
where
    R: AsyncRead + Unpin,
{
    let compression_method = read_compression_method(reader).await?;

    let content_type = read_content_type(reader).await?;
    validate_content_type(content_type)?;

    let _content_type_id = read_signed_int(reader, version).await?;
    let compressed_size: u64 = read_unsigned_int_as(reader, version).await?;
    let uncompressed_size: u64 = read_unsigned_int_as(reader, version).await?;

    let decoder = match compression_method {
        CompressionMethod::None => Decoder::None {
            inner: reader.take(uncompressed_size),
        },
        CompressionMethod::Gzip => Decoder::Gzip {
            inner: GzipDecoder::new(BufReader::new(reader.take(compressed_size))),
        },
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid block compression method: expected {:?} or {:?}, got {:?}",
                    CompressionMethod::None,
                    CompressionMethod::Gzip,
                    compression_method,
                ),
            ));
        }
    };

    // CRC32 bytes (if present for version >= 3.0) are consumed by the container's Take limit.

    Ok(decoder)
}

pin_project! {
    /// An async CRAM header container block data decoder.
    #[allow(missing_docs)]
    #[project = DecoderProjection]
    pub enum Decoder<R> {
        /// Uncompressed.
        None { #[pin] inner: Take<R> },
        /// gzip-compressed.
        Gzip { #[pin] inner: GzipDecoder<BufReader<Take<R>>> },
    }
}

impl<R> AsyncRead for Decoder<R>
where
    R: AsyncRead + Unpin,
{
    fn poll_read(
        self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &mut ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        match self.project() {
            DecoderProjection::None { inner } => inner.poll_read(cx, buf),
            DecoderProjection::Gzip { inner } => inner.poll_read(cx, buf),
        }
    }
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
