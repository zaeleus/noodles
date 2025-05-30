use std::io::{self, Read};

use flate2::Crc;

use crate::{BGZF_HEADER_SIZE, gz, io::Block};

const MIN_FRAME_SIZE: usize = BGZF_HEADER_SIZE + gz::TRAILER_SIZE;

type HeaderBuf = [u8; BGZF_HEADER_SIZE];
type TrailerBuf = [u8; gz::TRAILER_SIZE];

pub(crate) fn read_frame_into<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<Option<()>>
where
    R: Read,
{
    buf.resize(BGZF_HEADER_SIZE, 0);

    match reader.read_exact(buf) {
        Ok(()) => {}
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }

    // SAFETY: `buf.len() == BGZF_HEADER_SIZE >= mem::size_of::<u16>()`.
    let bsize = buf.last_chunk().map(|b| u16::from_le_bytes(*b)).unwrap();
    let block_size = usize::from(bsize) + 1;

    if block_size < MIN_FRAME_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid frame size",
        ));
    }

    buf.resize(block_size, 0);
    reader.read_exact(&mut buf[BGZF_HEADER_SIZE..])?;

    Ok(Some(()))
}

fn split_frame(buf: &[u8]) -> io::Result<(&HeaderBuf, &[u8], &TrailerBuf)> {
    if buf.len() < MIN_FRAME_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "invalid frame size",
        ));
    }

    // SAFETY: `buf.len() >= BGZF_HEADER_SIZE`.
    let (header, _) = buf.split_first_chunk().unwrap();

    let end = buf.len() - gz::TRAILER_SIZE;
    let cdata = &buf[BGZF_HEADER_SIZE..end];

    // SAFETY: `buf.len() >= gz::TRAILER_SIZE`.
    let (_, trailer) = buf.split_last_chunk().unwrap();

    Ok((header, cdata, trailer))
}

fn parse_header(src: &HeaderBuf) -> io::Result<()> {
    if is_valid_header(src) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ))
    }
}

fn is_valid_header(src: &HeaderBuf) -> bool {
    const BGZF_CM: u8 = 0x08; // DEFLATE
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XLEN: [u8; 2] = [0x06, 0x00];
    const BGZF_SI: [u8; 2] = [b'B', b'C'];
    const BGZF_SLEN: [u8; 2] = [0x02, 0x00];

    src[0..2] == gz::MAGIC_NUMBER
        && src[2] == BGZF_CM
        && src[3] == BGZF_FLG
        && src[10..12] == BGZF_XLEN
        && src[12..14] == BGZF_SI
        && src[14..16] == BGZF_SLEN
}

fn parse_trailer(src: &TrailerBuf) -> io::Result<(u32, usize)> {
    // SAFETY: `src.len() == 8`.
    let crc32 = u32::from_le_bytes(src[..4].try_into().unwrap());

    // SAFETY: `src.len() == 8`.
    let isize = usize::try_from(u32::from_le_bytes(src[4..].try_into().unwrap()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((crc32, isize))
}

pub(crate) fn parse_block(src: &[u8], block: &mut Block) -> io::Result<()> {
    let (block_size, cdata, crc32, isize) = parse_frame(src)?;
    block_initialize(block, block_size, isize);
    inflate(cdata, crc32, block.data_mut().as_mut())?;
    Ok(())
}

pub(super) fn parse_block_into_buf(
    src: &[u8],
    block: &mut Block,
    buf: &mut [u8],
) -> io::Result<()> {
    let (block_size, cdata, crc32, isize) = parse_frame(src)?;
    block_initialize(block, block_size, isize);
    block.data_mut().set_position(isize);
    inflate(cdata, crc32, &mut buf[..isize])?;
    Ok(())
}

fn parse_frame(src: &[u8]) -> io::Result<(u64, &[u8], u32, usize)> {
    let (header, cdata, trailer) = split_frame(src)?;

    let block_size =
        u64::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    parse_header(header)?;
    let (crc32, isize) = parse_trailer(trailer)?;

    Ok((block_size, cdata, crc32, isize))
}

fn block_initialize(block: &mut Block, block_size: u64, isize: usize) {
    block.set_size(block_size);

    let data = block.data_mut();
    data.set_position(0);
    data.resize(isize);
}

fn inflate(src: &[u8], crc32: u32, dst: &mut [u8]) -> io::Result<()> {
    use crate::deflate;

    deflate::decode(src, dst)?;

    let mut crc = Crc::new();
    crc.update(dst);

    if crc.sum() == crc32 {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data checksum mismatch",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::writer::BGZF_EOF;

    #[test]
    fn test_parse_header() -> io::Result<()> {
        parse_header(BGZF_EOF[0..BGZF_HEADER_SIZE].try_into().unwrap())?;
        Ok(())
    }

    #[test]
    fn test_is_valid_header() {
        let mut src = [
            0x1f, 0x8b, // ID1, ID2
            0x08, // CM = DEFLATE
            0x04, // FLG = FEXTRA
            0x00, 0x00, 0x00, 0x00, // MTIME = 0
            0x00, // XFL = 0
            0xff, // OS = 255 (unknown)
            0x06, 0x00, // XLEN = 6
            b'B', b'C', // SI1, SI2
            0x02, 0x00, // SLEN = 2
            0x1b, 0x00, // BSIZE = 27
        ];

        assert!(is_valid_header(&src));

        src[0] = 0x00;
        assert!(!is_valid_header(&src));
    }

    #[test]
    fn test_parse_trailer() -> io::Result<()> {
        let (_, src) = BGZF_EOF.split_last_chunk().unwrap();

        let (crc32, isize) = parse_trailer(src)?;
        assert_eq!(crc32, 0);
        assert_eq!(isize, 0);

        Ok(())
    }

    #[test]
    fn test_read_frame_into() -> Result<(), Box<dyn std::error::Error>> {
        let mut src = &BGZF_EOF[..];
        let mut buf = Vec::new();
        read_frame_into(&mut src, &mut buf)?.ok_or("invalid frame")?;
        assert_eq!(buf, BGZF_EOF);
        Ok(())
    }

    #[test]
    fn test_read_frame_into_with_invalid_block_size() {
        let data = {
            let mut eof = BGZF_EOF.to_vec();
            // BSIZE = 0
            eof[16] = 0x00;
            eof[17] = 0x00;
            eof
        };

        let mut reader = &data[..];
        let mut buf = Vec::new();
        assert!(read_frame_into(&mut reader, &mut buf).is_err());
    }
}
