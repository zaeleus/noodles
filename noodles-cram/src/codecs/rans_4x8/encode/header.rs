use std::io;

use crate::codecs::rans_4x8::Order;

pub const SIZE: usize = 9;

pub fn write_header(
    dst: &mut [u8; SIZE],
    order: Order,
    compressed_size: usize,
    uncompressed_size: usize,
) -> io::Result<()> {
    write_order(&mut dst[0], order);
    write_size(&mut dst[1..5], compressed_size)?;
    write_size(&mut dst[5..SIZE], uncompressed_size)?;
    Ok(())
}

fn write_order(dst: &mut u8, order: Order) {
    fn encode(order: Order) -> u8 {
        match order {
            Order::Zero => 0,
            Order::One => 1,
        }
    }

    *dst = encode(order);
}

fn write_size(dst: &mut [u8], size: usize) -> io::Result<()> {
    let n = u32::try_from(size).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.copy_from_slice(&n.to_le_bytes());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut dst = [0; SIZE];
        write_header(&mut dst, Order::Zero, 8, 13)?;

        let expected = [
            0x00, // order 0
            0x08, 0x00, 0x00, 0x00, // compressed size = 8
            0x0d, 0x00, 0x00, 0x00, // uncompressed size = 13
        ];

        assert_eq!(dst, expected);

        Ok(())
    }

    #[test]
    fn test_write_order() {
        let mut dst = 0;

        write_order(&mut dst, Order::Zero);
        assert_eq!(dst, 0);

        write_order(&mut dst, Order::One);
        assert_eq!(dst, 1);
    }
}
