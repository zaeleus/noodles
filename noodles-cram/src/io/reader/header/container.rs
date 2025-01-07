mod block;
mod header;

use std::io::{self, BufRead, BufReader, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use bytes::BytesMut;
use noodles_sam as sam;

use self::{block::read_block, header::read_header};

pub fn read_header_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<sam::Header>
where
    R: Read,
{
    let len = read_header(reader)?;

    buf.resize(len, 0);
    reader.read_exact(buf)?;

    let buf = buf.split().freeze();
    let mut reader = &buf[..];
    read_sam_header(&mut reader)
}

fn read_sam_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    let mut block_reader = read_block(reader)?;

    let len = block_reader.read_i32::<LittleEndian>().and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut parser = sam::header::Parser::default();

    let mut header_reader = BufReader::new(block_reader.take(len));
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
    use super::*;

    #[test]
    fn test_read_sam_header() -> io::Result<()> {
        use sam::header::record::value::{
            map::{self, header::Version},
            Map,
        };

        let mut src = vec![
            0x00, // compression method = none (0)
            0x00, // content type = file header (0)
            0x00, // block content ID
            0x0f, // compressed size = 15
            0x0f, // uncompressed size = 15
        ];
        src.extend(11i32.to_le_bytes());
        src.extend(b"@HD\tVN:1.6\n");
        src.extend([0x00, 0x00, 0x00, 0x00]); // CRC32

        let mut reader = &src[..];
        let actual = read_sam_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_sam_header_with_invalid_compression_method() {
        let src = [
            0x03, // compression method = LZMA (3)
            0x00, // content type = file header (0)
            0x00, // block content ID
            0x0f, // compressed size = 15
            0x0f, // uncompressed size = 15
            // ...
            0x00, 0x00, 0x00, 0x00, // CRC32
        ];

        let mut reader = &src[..];
        assert!(matches!(
            read_sam_header(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_read_sam_header_with_invalid_content_type() {
        let src = [
            0x00, // compression method = none (0)
            0x04, // content type = external data (4)
        ];

        let mut reader = &src[..];

        assert!(matches!(
            read_sam_header(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
