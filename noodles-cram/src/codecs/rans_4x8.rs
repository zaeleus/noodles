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
    fn test_self() -> io::Result<()> {
        let src = b"noodles";
        let options = [Order::Zero, Order::One];

        for order in options {
            let compressed_data = encode(order, src)?;
            let uncompressed_data = decode(&mut &compressed_data[..])?;
            assert_eq!(uncompressed_data, src);
        }

        Ok(())
    }
}
