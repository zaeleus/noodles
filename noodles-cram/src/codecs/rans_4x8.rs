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
    use crate::codecs::test_utils::sparse;

    #[test]
    fn test_decode_htscodecs_fixtures() -> io::Result<()> {
        let qual = include_bytes!("fixtures/qual.bin");
        assert_eq!(decode(include_bytes!("fixtures/r4x8-qual-o0.bin"))?, qual);
        assert_eq!(decode(include_bytes!("fixtures/r4x8-qual-o1.bin"))?, qual);
        assert!(decode(include_bytes!("fixtures/r4x8-big-o1.bin"))? == sparse(70001, 2048, 4));
        assert!(decode(include_bytes!("fixtures/r4x8-wide-o0.bin"))? == sparse(4099, 16, 32));
        Ok(())
    }

    #[test]
    fn test_self() -> io::Result<()> {
        let src = b"noodles";
        let options = [Order::Zero, Order::One];

        for order in options {
            let compressed_data = encode(order, src)?;
            let uncompressed_data = decode(&compressed_data)?;
            assert_eq!(uncompressed_data, src);
        }

        Ok(())
    }
}
