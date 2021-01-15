mod decode;
mod encode;
mod order;

pub use self::{decode::rans_decode, encode::rans_encode, order::Order};

#[cfg(test)]
mod tests {
    use std::io;

    use super::*;

    #[test]
    fn test_self_0() -> io::Result<()> {
        let data = b"noodles";

        let compressed_data = rans_encode(Order::Zero, data)?;

        let mut reader = &compressed_data[..];
        let decompressed_data = rans_decode(&mut reader)?;

        assert_eq!(decompressed_data, data);

        Ok(())
    }

    #[test]
    fn test_self_1() -> io::Result<()> {
        let data = b"abracadabraabracadabraabracadabraabracadabra";

        let compressed_data = rans_encode(Order::One, data)?;

        let mut reader = &compressed_data[..];
        let decompressed_data = rans_decode(&mut reader)?;

        assert_eq!(decompressed_data, data);

        Ok(())
    }
}
