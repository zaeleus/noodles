use std::io::{self, Read};

use byteorder::{ByteOrder, LittleEndian, ReadBytesExt};
use flate2::Crc;

use super::inflate_data;
use crate::{gz, Block, BGZF_HEADER_SIZE};

fn read_header<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut header = [0; BGZF_HEADER_SIZE];

    match reader.read_exact(&mut header) {
        Ok(_) => {}
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    }

    if !is_valid_header(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ));
    }

    let bsize = LittleEndian::read_u16(&header[16..]);

    // Add 1 because BSIZE is "total Block SIZE minus 1".
    Ok(u32::from(bsize) + 1)
}

fn is_valid_header(header: &[u8; BGZF_HEADER_SIZE]) -> bool {
    const BGZF_CM: u8 = 0x08; // DEFLATE
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XLEN: u16 = 6;
    const BGZF_SI1: u8 = b'B';
    const BGZF_SI2: u8 = b'C';
    const BGZF_SLEN: u16 = 2;

    let magic_number = &header[0..2];
    let cm = header[2];
    let flg = header[3];
    let xlen = LittleEndian::read_u16(&header[10..]);
    let subfield_id_1 = header[12];
    let subfield_id_2 = header[13];
    let bc_len = LittleEndian::read_u16(&header[14..]);

    magic_number == gz::MAGIC_NUMBER
        && cm == BGZF_CM
        && flg == BGZF_FLG
        && xlen == BGZF_XLEN
        && subfield_id_1 == BGZF_SI1
        && subfield_id_2 == BGZF_SI2
        && bc_len == BGZF_SLEN
}

/// Reads a BGZF block trailer.
///
/// The position of the stream is expected to be at the start of the block trailer, i.e., 8 bytes
/// from the end of the block.
fn read_trailer<R>(reader: &mut R) -> io::Result<(u32, usize)>
where
    R: Read,
{
    let crc32 = reader.read_u32::<LittleEndian>()?;

    let r#isize = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    Ok((crc32, r#isize))
}

fn read_compressed_block<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<(usize, (u32, usize))>
where
    R: Read,
{
    let clen = match read_header(reader) {
        Ok(0) => return Ok((0, (0, 0))),
        Ok(bs) => bs as usize,
        Err(e) => return Err(e),
    };

    if clen < BGZF_HEADER_SIZE + gz::TRAILER_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "expected clen >= {}, got {}",
                BGZF_HEADER_SIZE + gz::TRAILER_SIZE,
                clen
            ),
        ));
    }

    let cdata_len = clen - BGZF_HEADER_SIZE - gz::TRAILER_SIZE;
    buf.resize(cdata_len, Default::default());
    reader.read_exact(buf)?;

    let trailer = read_trailer(reader)?;

    Ok((clen, trailer))
}

pub(super) fn read_block<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    block: &mut Block,
) -> io::Result<usize>
where
    R: Read,
{
    let (clen, crc32, ulen) = match read_compressed_block(reader, buf) {
        Ok((0, (_, 0))) => return Ok(0),
        Ok((clen, (crc32, ulen))) => (clen, crc32, ulen),
        Err(e) => return Err(e),
    };

    block.set_size(clen as u64);

    let data = block.data_mut();
    data.set_position(0);
    data.resize(ulen);

    inflate_data(buf, data.as_mut())?;

    let mut crc = Crc::new();
    crc.update(data.as_ref());

    if crc.sum() == crc32 {
        Ok(clen)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data checksum mismatch",
        ))
    }
}

pub(super) fn read_block_into<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    block: &mut Block,
    dst: &mut [u8],
) -> io::Result<usize>
where
    R: Read,
{
    let (clen, crc32, ulen) = match read_compressed_block(reader, buf) {
        Ok((0, (_, 0))) => return Ok(0),
        Ok((clen, (crc32, ulen))) => (clen, crc32, ulen),
        Err(e) => return Err(e),
    };

    block.set_size(clen as u64);

    let data = block.data_mut();
    data.resize(ulen);
    data.set_position(ulen);

    inflate_data(buf, &mut dst[..ulen])?;

    let mut crc = Crc::new();
    crc.update(&dst[..ulen]);

    if crc.sum() == crc32 {
        Ok(clen)
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
    fn test_read_header() -> io::Result<()> {
        let mut reader = BGZF_EOF;
        let block_size = read_header(&mut reader)?;
        assert_eq!(block_size, BGZF_EOF.len() as u32);
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
    fn test_read_trailer() -> io::Result<()> {
        let (_, mut reader) = BGZF_EOF.split_at(BGZF_EOF.len() - gz::TRAILER_SIZE);

        let (crc32, r#isize) = read_trailer(&mut reader)?;
        assert_eq!(crc32, 0);
        assert_eq!(r#isize, 0);

        Ok(())
    }

    #[test]
    fn test_read_trailer_with_invalid_block_trailer() {
        let data = [0, 0, 0, 0];
        let mut reader = &data[..];
        assert!(read_trailer(&mut reader).is_err());
    }

    #[test]
    fn test_read_block() -> io::Result<()> {
        let mut reader = BGZF_EOF;
        let mut cdata = Vec::new();
        let mut block = Block::default();

        let block_size = read_block(&mut reader, &mut cdata, &mut block)?;
        assert_eq!(block_size, BGZF_EOF.len());

        Ok(())
    }

    #[test]
    fn test_read_block_with_invalid_block_size() {
        let data = {
            let mut eof = BGZF_EOF.to_vec();
            // BSIZE = 0
            eof[16] = 0x00;
            eof[17] = 0x00;
            eof
        };

        let mut reader = &data[..];
        let mut cdata = Vec::new();
        let mut block = Block::default();

        assert!(read_block(&mut reader, &mut cdata, &mut block).is_err());
    }
}
