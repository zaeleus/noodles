pub mod multi;
pub mod single;

use std::io::{self, Read};

use bytes::Buf;
use flate2::Crc;

use crate::{gz, Block, BGZF_HEADER_SIZE};

pub enum Inner<R> {
    Single(single::Reader<R>),
    Multi(multi::Reader<R>),
}

impl<R> Inner<R>
where
    R: Read,
{
    pub fn get_ref(&self) -> &R {
        match self {
            Self::Single(reader) => reader.get_ref(),
            Self::Multi(reader) => reader.get_ref(),
        }
    }

    pub fn get_mut(&mut self) -> &mut R {
        match self {
            Self::Single(reader) => reader.get_mut(),
            Self::Multi(reader) => reader.get_mut(),
        }
    }

    pub fn into_inner(self) -> R {
        match self {
            Self::Single(reader) => reader.into_inner(),
            Self::Multi(reader) => reader.into_inner(),
        }
    }

    pub fn next_block(&mut self) -> io::Result<Option<Block>> {
        match self {
            Self::Single(reader) => reader.next_block(),
            Self::Multi(reader) => reader.next_block(),
        }
    }
}

fn read_frame<R>(reader: &mut R) -> io::Result<Option<Vec<u8>>>
where
    R: Read,
{
    let mut buf = vec![0; BGZF_HEADER_SIZE];

    if read_frame_into(reader, &mut buf)?.is_some() {
        Ok(Some(buf))
    } else {
        Ok(None)
    }
}

fn read_frame_into<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<Option<()>>
where
    R: Read,
{
    const MIN_FRAME_SIZE: usize = BGZF_HEADER_SIZE + gz::TRAILER_SIZE;
    const BSIZE_POSITION: usize = 16;

    buf.resize(BGZF_HEADER_SIZE, 0);

    match reader.read_exact(buf) {
        Ok(()) => {}
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }

    let bsize = (&buf[BSIZE_POSITION..]).get_u16_le();
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

fn split_frame(buf: &[u8]) -> (&[u8], &[u8], &[u8]) {
    let header = &buf[..BGZF_HEADER_SIZE];

    let n = buf.len() - gz::TRAILER_SIZE;
    let cdata = &buf[BGZF_HEADER_SIZE..n];

    let trailer = &buf[n..];

    (header, cdata, trailer)
}

fn parse_header(src: &[u8]) -> io::Result<()> {
    if is_valid_header(src) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ))
    }
}

fn is_valid_header<B>(mut src: B) -> bool
where
    B: Buf,
{
    use std::mem;

    const BGZF_CM: u8 = 0x08; // DEFLATE
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XLEN: u16 = 6;
    const BGZF_SI1: u8 = b'B';
    const BGZF_SI2: u8 = b'C';
    const BGZF_SLEN: u16 = 2;

    let id_1 = src.get_u8();
    let id_2 = src.get_u8();
    let cm = src.get_u8();
    let flg = src.get_u8();

    // 4 (MTIME) + 1 (XFL) + 1 (OS)
    src.advance(mem::size_of::<u32>() + mem::size_of::<u8>() + mem::size_of::<u8>());

    let xlen = src.get_u16_le();
    let subfield_id_1 = src.get_u8();
    let subfield_id_2 = src.get_u8();
    let subfield_len = src.get_u16_le();

    id_1 == gz::MAGIC_NUMBER[0]
        && id_2 == gz::MAGIC_NUMBER[1]
        && cm == BGZF_CM
        && flg == BGZF_FLG
        && xlen == BGZF_XLEN
        && subfield_id_1 == BGZF_SI1
        && subfield_id_2 == BGZF_SI2
        && subfield_len == BGZF_SLEN
}

fn parse_trailer<B>(mut src: B) -> io::Result<(u32, usize)>
where
    B: Buf,
{
    let crc32 = src.get_u32_le();

    let r#isize = usize::try_from(src.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((crc32, r#isize))
}

pub(crate) fn parse_frame(src: &[u8]) -> io::Result<Block> {
    let (header, cdata, trailer) = split_frame(src);

    parse_header(header)?;
    let (crc32, r#isize) = parse_trailer(trailer)?;

    let mut block = Block::default();

    let block_size =
        u64::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    block.set_size(block_size);

    let data = block.data_mut();
    data.set_position(0);
    data.resize(r#isize);

    inflate(cdata, crc32, data.as_mut())?;

    Ok(block)
}

fn inflate(src: &[u8], crc32: u32, dst: &mut [u8]) -> io::Result<()> {
    use super::inflate_data;

    inflate_data(src, dst)?;

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
    use crate::writer::BGZF_EOF;

    #[test]
    fn test_parse_header() -> io::Result<()> {
        parse_header(BGZF_EOF)?;
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

        let mut reader = &src[..];
        assert!(is_valid_header(&mut reader));

        src[0] = 0x00;
        let mut reader = &src[..];
        assert!(!is_valid_header(&mut reader));
    }

    #[test]
    fn test_parse_trailer() -> io::Result<()> {
        let (_, mut src) = BGZF_EOF.split_at(BGZF_EOF.len() - gz::TRAILER_SIZE);

        let (crc32, r#isize) = parse_trailer(&mut src)?;
        assert_eq!(crc32, 0);
        assert_eq!(r#isize, 0);

        Ok(())
    }

    #[test]
    fn test_read_frame() -> Result<(), Box<dyn std::error::Error>> {
        let mut src = BGZF_EOF;
        let buf = read_frame(&mut src)?.ok_or("invalid frame")?;
        assert_eq!(buf, BGZF_EOF);
        Ok(())
    }

    #[test]
    fn test_read_frame_with_invalid_block_size() {
        let data = {
            let mut eof = BGZF_EOF.to_vec();
            // BSIZE = 0
            eof[16] = 0x00;
            eof[17] = 0x00;
            eof
        };

        let mut reader = &data[..];
        assert!(read_frame(&mut reader).is_err());
    }
}
