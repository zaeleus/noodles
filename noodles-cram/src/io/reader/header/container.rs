mod header;

use std::io::{self, BufRead, BufReader, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use bytes::{Bytes, BytesMut};
use noodles_sam as sam;

use self::header::read_header;
use crate::container::{
    block::{CompressionMethod, ContentType},
    Block,
};

pub fn read_header_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<sam::Header>
where
    R: Read,
{
    let len = read_header(reader)?;

    buf.resize(len, 0);
    reader.read_exact(buf)?;
    let mut buf = buf.split().freeze();

    read_sam_header_from_block(&mut buf)
}

pub fn read_sam_header_from_block(src: &mut Bytes) -> io::Result<sam::Header> {
    use crate::io::reader::container::read_block;

    let block = read_block(src)?;
    read_sam_header(&block)
}

fn read_sam_header(block: &Block) -> io::Result<sam::Header> {
    use flate2::bufread::GzDecoder;

    const EXPECTED_CONTENT_TYPE: ContentType = ContentType::FileHeader;

    if block.content_type() != EXPECTED_CONTENT_TYPE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block content type: expected {:?}, got {:?}",
                EXPECTED_CONTENT_TYPE,
                block.content_type()
            ),
        ));
    }

    let mut reader: Box<dyn BufRead> = match block.compression_method() {
        CompressionMethod::None => Box::new(block.data()),
        CompressionMethod::Gzip => {
            let decoder = GzDecoder::new(block.data());
            Box::new(BufReader::new(decoder))
        }
        compression_method => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid block compression method: expected {:?} or {:?}, got {:?}",
                    CompressionMethod::None,
                    CompressionMethod::Gzip,
                    compression_method
                ),
            ))
        }
    };

    let len = reader.read_i32::<LittleEndian>().and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut parser = sam::header::Parser::default();

    let mut header_reader = reader.take(len);
    let mut buf = Vec::new();

    while read_line(&mut header_reader, &mut buf)? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    Ok(parser.finish())
}

fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst)? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use bytes::BufMut;

    use super::*;

    #[test]
    fn test_read_sam_header() -> io::Result<()> {
        use sam::header::record::value::{
            map::{self, header::Version},
            Map,
        };

        let raw_header = "@HD\tVN:1.6\n@CO\tnoodles-cram\n";

        let header_data = raw_header.to_string().into_bytes();
        let header_data_len = i32::try_from(header_data.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let mut data = Vec::new();
        data.put_i32_le(header_data_len);
        data.extend(&header_data);

        let block = Block::builder()
            .set_content_type(ContentType::FileHeader)
            .set_uncompressed_len(data.len())
            .set_data(data.into())
            .build();

        let actual = read_sam_header(&block)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_comment("noodles-cram")
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_sam_header_with_invalid_compression_method() {
        let block = Block::builder()
            .set_compression_method(CompressionMethod::Lzma)
            .set_content_type(ContentType::FileHeader)
            .build();

        assert!(matches!(
            read_sam_header(&block),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_read_sam_header_with_invalid_content_type() {
        let block = Block::builder()
            .set_content_type(ContentType::ExternalData)
            .build();

        assert!(matches!(
            read_sam_header(&block),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
