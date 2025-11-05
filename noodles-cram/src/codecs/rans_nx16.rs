//! rANS Nx16 codec.

pub(crate) mod decode;
pub(crate) mod encode;
mod flags;

pub use self::flags::Flags;
pub(crate) use self::{decode::decode, encode::encode};

const ALPHABET_SIZE: usize = 256;

#[cfg(test)]
mod tests {
    use std::io;

    use super::*;

    #[test]
    fn test_self() -> io::Result<()> {
        let src = b"noodles";

        let options = [
            Flags::empty(),
            Flags::ORDER,
            Flags::N32,
            Flags::STRIPE,
            Flags::CAT,
            // Flags::RLE,
            Flags::PACK,
        ];

        for flags in options {
            let compressed_data = encode(flags, src)?;
            let uncompressed_data = decode(&compressed_data, src.len())?;
            assert_eq!(uncompressed_data, src);
        }

        Ok(())
    }
}
