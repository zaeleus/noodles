use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::codecs::rans_4x8::Order;

pub fn write_header<W>(
    writer: &mut W,
    order: Order,
    compressed_size: usize,
    uncompressed_size: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_order(writer, order)?;
    write_size(writer, compressed_size)?;
    write_size(writer, uncompressed_size)?;
    Ok(())
}

fn write_order<W>(writer: &mut W, order: Order) -> io::Result<()>
where
    W: Write,
{
    writer.write_u8(u8::from(order))
}

fn write_size<W>(writer: &mut W, size: usize) -> io::Result<()>
where
    W: Write,
{
    let n = u32::try_from(size).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(n)
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
