use std::io::{self, Write};

use crate::{BGZF_HEADER_SIZE, gz};

pub(crate) fn write_frame<W>(
    writer: &mut W,
    compressed_data: &[u8],
    crc32: u32,
    uncompressed_size: usize,
) -> io::Result<usize>
where
    W: Write,
{
    let block_size = BGZF_HEADER_SIZE + compressed_data.len() + gz::TRAILER_SIZE;
    write_header(writer, block_size)?;

    writer.write_all(compressed_data)?;

    write_trailer(writer, crc32, uncompressed_size)?;

    Ok(block_size)
}

fn write_header<W>(writer: &mut W, block_size: usize) -> io::Result<()>
where
    W: Write,
{
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XFL: u8 = 0x00; // none
    const BGZF_XLEN: u16 = 6;

    const BGZF_SI1: u8 = b'B';
    const BGZF_SI2: u8 = b'C';
    const BGZF_SLEN: u16 = 2;

    writer.write_all(&gz::MAGIC_NUMBER)?;
    write_u8(writer, gz::CompressionMethod::Deflate as u8)?;
    write_u8(writer, BGZF_FLG)?;
    write_u32_le(writer, gz::MTIME_NONE)?;
    write_u8(writer, BGZF_XFL)?;
    write_u8(writer, gz::OperatingSystem::Unknown as u8)?;
    write_u16_le(writer, BGZF_XLEN)?;

    write_u8(writer, BGZF_SI1)?;
    write_u8(writer, BGZF_SI2)?;
    write_u16_le(writer, BGZF_SLEN)?;

    let bsize = u16::try_from(block_size - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u16_le(writer, bsize)?;

    Ok(())
}

fn write_trailer<W>(writer: &mut W, checksum: u32, uncompressed_size: usize) -> io::Result<()>
where
    W: Write,
{
    write_u32_le(writer, checksum)?;

    let isize = u32::try_from(uncompressed_size)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, isize)?;

    Ok(())
}

fn write_u8<W>(writer: &mut W, n: u8) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[n])
}

fn write_u16_le<W>(writer: &mut W, n: u16) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}

fn write_u32_le<W>(writer: &mut W, n: u32) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut buf = Vec::new();
        write_header(&mut buf, 8)?;

        let expected = [
            0x1f, 0x8b, // magic number
            0x08, // compression method = DEFLATE
            0x04, // flags = FEXTRA
            0x00, 0x00, 0x00, 0x00, // modification time = 0
            0x00, // extra flags = none
            0xff, // operating system = unknown
            0x06, 0x00, // extra subfields size = 6
            b'B', b'C', // subfields[0].id = b"BC"
            0x02, 0x00, // subfields[0].size = 2
            0x07, 0x00, // subfields[0].block_size = 7
        ];

        assert_eq!(buf, expected);

        assert!(matches!(
            write_header(&mut io::sink(), (1 << 16) + 1),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_trailer() -> io::Result<()> {
        let mut buf = Vec::new();

        write_trailer(&mut buf, 0x05080d15, 34)?;

        let expected = [
            0x15, 0x0d, 0x08, 0x05, // CRC32 = 0x05080d15
            0x22, 0x00, 0x00, 0x00, // ISIZE = 34
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
