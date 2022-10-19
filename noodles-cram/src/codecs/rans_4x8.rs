//! rANS 4x8 codec.

mod decode;
mod encode;
mod order;

pub use self::order::Order;
pub(crate) use self::{decode::decode, encode::encode};

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
