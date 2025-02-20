//! rANS 4x8 codec.

mod decode;
mod encode;
mod order;

pub use self::order::Order;
pub(crate) use self::{decode::decode, encode::encode};

// § 2.2 "rANS entropy encoding" (2023-03-15)
const ALPHABET_SIZE: usize = 256; // b
const STATE_COUNT: usize = 4; // |R|
const LOWER_BOUND: u32 = 1 << 23; // L

#[cfg(test)]
mod tests {
    use std::io;

    use super::*;

    #[test]
    fn test_self_0() -> io::Result<()> {
        let data = b"noodles";

        let compressed_data = encode(Order::Zero, data)?;

        let mut reader = &compressed_data[..];
        let decompressed_data = decode(&mut reader)?;

        assert_eq!(decompressed_data, data);

        Ok(())
    }

    #[test]
    fn test_self_1() -> io::Result<()> {
        let data = b"noodles";

        let compressed_data = encode(Order::One, data)?;

        let mut reader = &compressed_data[..];
        let decompressed_data = decode(&mut reader)?;

        assert_eq!(decompressed_data, data);

        Ok(())
    }
}
