mod header;
mod order_0;
mod order_1;

use std::io;

use self::header::read_header;
use super::{LOWER_BOUND, Order, STATE_COUNT};
use crate::io::reader::num::{read_u8, read_u32_le};

pub fn decode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let (order, _, uncompressed_size) = read_header(&mut src)?;

    let mut dst = vec![0; uncompressed_size];

    match order {
        Order::Zero => order_0::decode(&mut src, &mut dst)?,
        Order::One => order_1::decode(&mut src, &mut dst)?,
    }

    Ok(dst)
}

fn read_states(src: &mut &[u8]) -> io::Result<[u32; 4]> {
    let mut states = [0; STATE_COUNT];

    for state in &mut states {
        *state = read_u32_le(src)?;
    }

    Ok(states)
}

pub fn state_cumulative_frequency(s: u32) -> u16 {
    (s & 0x0fff) as u16
}

pub fn state_step(s: u32, f: u16, g: u16) -> u32 {
    u32::from(f) * (s >> 12) + (s & 0x0fff) - u32::from(g)
}

pub fn state_renormalize(mut s: u32, src: &mut &[u8]) -> io::Result<u32> {
    while s < LOWER_BOUND {
        let b = read_u8(src).map(u32::from)?;
        s = (s << 8) | b;
    }

    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_header() -> io::Result<()> {
        let src = [
            0x00, // order = 0
            0x25, 0x00, 0x00, 0x00, // compressed size = 37
            0x07, 0x00, 0x00, 0x00, // data size = 7
        ];

        assert_eq!(read_header(&mut &src[..])?, (Order::Zero, 37, 7));

        Ok(())
    }

    #[test]
    fn test_decode_with_order_0() -> io::Result<()> {
        let expected = b"noodles";

        let src = [
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49, 0x65, 0x00,
            0x82, 0x49, 0x6c, 0x82, 0x49, 0x6e, 0x82, 0x49, 0x6f, 0x00, 0x84, 0x92, 0x73, 0x82,
            0x49, 0x00, 0xe2, 0x06, 0x83, 0x18, 0x74, 0x7b, 0x41, 0x0c, 0x2b, 0xa9, 0x41, 0x0c,
            0x25, 0x31, 0x80, 0x03,
        ];

        let actual = decode(&src)?;
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_decode_with_order_1() -> io::Result<()> {
        let expected = b"noodles";

        let src = [
            0x01, 0x3b, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x64, 0x84, 0x00, 0x6e,
            0x84, 0x00, 0x6f, 0x00, 0x87, 0xff, 0x00, 0x64, 0x6c, 0x8f, 0xff, 0x00, 0x65, 0x00,
            0x73, 0x8f, 0xff, 0x00, 0x6c, 0x65, 0x8f, 0xff, 0x00, 0x6e, 0x6f, 0x8f, 0xff, 0x00,
            0x6f, 0x00, 0x64, 0x87, 0xff, 0x6f, 0x88, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x02,
            0x02, 0x28, 0x00, 0x01, 0x02, 0x28, 0x00, 0x01, 0x02, 0x60, 0x00, 0x02,
        ];

        let actual = decode(&src)?;
        assert_eq!(actual, expected);

        Ok(())
    }
}
