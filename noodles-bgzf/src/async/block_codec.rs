use std::io;

use bytes::{Buf, BufMut, Bytes, BytesMut};
use tokio_util::codec::{Decoder, Encoder};

use crate::{gz, BGZF_HEADER_SIZE};

use super::writer::deflate::GzData;

pub struct BlockCodec;

impl Decoder for BlockCodec {
    type Item = Bytes;
    type Error = io::Error;

    fn decode(&mut self, src: &mut BytesMut) -> Result<Option<Self::Item>, Self::Error> {
        if src.len() < BGZF_HEADER_SIZE {
            src.reserve(BGZF_HEADER_SIZE);
            return Ok(None);
        }

        let block_size = {
            let mut header = &src[..BGZF_HEADER_SIZE];
            header.advance(16); // [ID1, ..., SLEN]
            usize::from(header.get_u16_le()) + 1
        };

        if src.len() < block_size {
            src.reserve(block_size);
            return Ok(None);
        }

        Ok(Some(src.split_to(block_size).freeze()))
    }
}

impl Encoder<GzData> for BlockCodec {
    type Error = io::Error;

    fn encode(
        &mut self,
        (cdata, crc32, r#isize): GzData,
        dst: &mut BytesMut,
    ) -> Result<(), Self::Error> {
        let cdata_len = cdata.len();
        let block_size = BGZF_HEADER_SIZE + cdata_len + gz::TRAILER_SIZE;

        dst.reserve(block_size);

        put_header(dst, block_size)?;
        dst.extend_from_slice(&cdata);
        put_trailer(dst, crc32, r#isize);

        Ok(())
    }
}

fn put_header(dst: &mut BytesMut, block_size: usize) -> io::Result<()> {
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XFL: u8 = 0x00; // none
    const BGZF_XLEN: u16 = 6;

    const BGZF_SI1: u8 = 0x42;
    const BGZF_SI2: u8 = 0x43;
    const BGZF_SLEN: u16 = 2;

    dst.extend_from_slice(&gz::MAGIC_NUMBER);
    dst.put_u8(gz::CompressionMethod::Deflate as u8);
    dst.put_u8(BGZF_FLG);
    dst.put_u32_le(gz::MTIME_NONE);
    dst.put_u8(BGZF_XFL);
    dst.put_u8(gz::OperatingSystem::Unknown as u8);
    dst.put_u16_le(BGZF_XLEN);

    dst.put_u8(BGZF_SI1);
    dst.put_u8(BGZF_SI2);
    dst.put_u16_le(BGZF_SLEN);

    let bsize = u16::try_from(block_size - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u16_le(bsize);

    Ok(())
}

fn put_trailer(dst: &mut BytesMut, checksum: u32, uncompressed_size: u32) {
    dst.put_u32_le(checksum);
    dst.put_u32_le(uncompressed_size);
}

#[cfg(test)]
mod tests {
    use super::*;

    static BLOCK: &[u8] = &[
        0x1f, 0x8b, // (ID1, ID2)
        0x08, // CM = DEFLATE
        0x04, // FLG = FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME = 0
        0x00, // XFL = 0
        0xff, // OS = 255 (unknown)
        0x06, 0x00, // XLEN = 6
        0x42, 0x43, // (SI1, SI2)
        0x02, 0x00, // SLEN = 2
        0x22, 0x00, // BSIZE = 34 (+ 1)
        0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, // CDATA = deflate(b"noodles")
        0xa1, 0x58, 0x2a, 0x80, // CRC32
        0x07, 0x00, 0x00, 0x00, // ISIZE
    ];

    #[test]
    fn test_decode() -> io::Result<()> {
        let mut decoder = BlockCodec;

        let mut src = BytesMut::from(BLOCK);

        let block = decoder.decode(&mut src)?;
        assert_eq!(block.as_deref(), Some(BLOCK));

        let block = decoder.decode(&mut src)?;
        assert!(block.is_none());

        Ok(())
    }

    #[test]
    fn test_encode() -> io::Result<()> {
        let mut encoder = BlockCodec;

        let cdata = vec![0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00];
        let crc32 = 0x802a58a1;
        let r#isize = 7;
        let data = (cdata, crc32, r#isize);

        let mut dst = BytesMut::new();

        encoder.encode(data, &mut dst)?;

        assert_eq!(&dst[..], BLOCK);

        Ok(())
    }

    #[test]
    fn test_put_header() {
        let mut dst = BytesMut::new();

        assert!(put_header(&mut dst, 8).is_ok());

        assert!(matches!(
            put_header(&mut dst, (1 << 16) + 1),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));
    }
}
