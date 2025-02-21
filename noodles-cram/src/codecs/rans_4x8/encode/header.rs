use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::codecs::rans_4x8::Order;

pub fn write_header<W>(
    writer: &mut W,
    order: Order,
    compressed_len: usize,
    uncompressed_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    writer.write_u8(u8::from(order))?;

    let compressed_size = u32::try_from(compressed_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(compressed_size)?;

    let data_size = u32::try_from(uncompressed_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(data_size)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut writer = Vec::new();
        write_header(&mut writer, Order::One, 14930352, 9227465)?;

        let expected = [
            0x01, // order
            0xb0, 0xd1, 0xe3, 0x00, // compressed length
            0xc9, 0xcc, 0x8c, 0x00, // data length
        ];

        assert_eq!(writer, expected);

        Ok(())
    }
}
