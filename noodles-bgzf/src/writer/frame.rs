use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::{gz, BGZF_HEADER_SIZE};

pub(crate) fn write_frame<W>(
    writer: &mut W,
    compressed_data: &[u8],
    crc32: u32,
    uncompressed_len: usize,
) -> io::Result<usize>
where
    W: Write,
{
    let block_size = BGZF_HEADER_SIZE + compressed_data.len() + gz::TRAILER_SIZE;
    write_header(writer, block_size)?;

    writer.write_all(compressed_data)?;

    write_trailer(writer, crc32, uncompressed_len)?;

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
    writer.write_u8(gz::CompressionMethod::Deflate as u8)?;
    writer.write_u8(BGZF_FLG)?;
    writer.write_u32::<LittleEndian>(gz::MTIME_NONE)?;
    writer.write_u8(BGZF_XFL)?;
    writer.write_u8(gz::OperatingSystem::Unknown as u8)?;
    writer.write_u16::<LittleEndian>(BGZF_XLEN)?;

    writer.write_u8(BGZF_SI1)?;
    writer.write_u8(BGZF_SI2)?;
    writer.write_u16::<LittleEndian>(BGZF_SLEN)?;

    let bsize = u16::try_from(block_size - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(bsize)?;

    Ok(())
}

fn write_trailer<W>(writer: &mut W, checksum: u32, uncompressed_len: usize) -> io::Result<()>
where
    W: Write,
{
    writer.write_u32::<LittleEndian>(checksum)?;

    let r#isize = u32::try_from(uncompressed_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(r#isize)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() {
        let mut writer = io::sink();

        assert!(write_header(&mut writer, 8).is_ok());

        assert!(matches!(
            write_header(&mut writer, (1 << 16) + 1),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));
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
