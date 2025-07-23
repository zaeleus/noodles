use std::io::{self, Write};

use crate::{
    codecs::rans_4x8::Order,
    io::writer::num::{write_u8, write_u32_le},
};

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
    let n = match order {
        Order::Zero => 0,
        Order::One => 1,
    };

    write_u8(writer, n)
}

fn write_size<W>(writer: &mut W, size: usize) -> io::Result<()>
where
    W: Write,
{
    let n = u32::try_from(size).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut buf = Vec::new();
        write_header(&mut buf, Order::Zero, 8, 13)?;

        let expected = [
            0x00, // order 0
            0x08, 0x00, 0x00, 0x00, // compressed size = 8
            0x0d, 0x00, 0x00, 0x00, // uncompressed size = 13
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_order() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_order(&mut buf, Order::Zero)?;
        assert_eq!(buf, [0x00]);

        buf.clear();
        write_order(&mut buf, Order::One)?;
        assert_eq!(buf, [0x01]);

        Ok(())
    }
}
