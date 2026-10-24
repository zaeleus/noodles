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
    use crate::codecs::test_utils::sparse;

    #[test]
    fn test_decode_htscodecs_fixtures() -> io::Result<()> {
        macro_rules! r4x16 {
            ($name:literal) => {
                (
                    $name,
                    &include_bytes!(concat!("fixtures/r4x16-", $name, ".bin"))[..],
                )
            };
        }

        let qual = include_bytes!("fixtures/qual.bin").to_vec();
        // `mid` has 12-bit order-1 frequencies.
        let (mid, big) = (sparse(8000, 1024, 8), sparse(70001, 2048, 4));

        // Named by flags: ORDER = 0x01, N32 = 0x04, STRIPE = 0x08, CAT = 0x20, RLE = 0x40 and
        // PACK = 0x80.
        for (name, src) in [
            r4x16!("qual-00"),
            r4x16!("qual-01"),
            r4x16!("qual-04"),
            r4x16!("qual-05"),
            r4x16!("qual-41"),
            r4x16!("qual-49"),
            r4x16!("qual-64"),
            r4x16!("qual-80"),
            r4x16!("qual-c0"),
            r4x16!("mid-01"),
            r4x16!("big-05"),
        ] {
            let expected = match &name[..3] {
                "qua" => &qual,
                "mid" => &mid,
                _ => &big,
            };

            assert_eq!(format!("{:02x}", src[0]), name[name.len() - 2..]);
            assert!(decode(src, expected.len())? == *expected, "{name}");
        }

        Ok(())
    }

    #[test]
    fn test_self() -> io::Result<()> {
        let src = b"noodles";

        let options = [
            Flags::empty(),
            Flags::ORDER,
            Flags::N32,
            Flags::STRIPE,
            Flags::CAT,
            Flags::RLE,
            Flags::PACK,
        ];

        for flags in options {
            let compressed_data = encode(flags, src)?;
            let uncompressed_data = decode(&compressed_data, src.len())?;
            assert_eq!(uncompressed_data, src);
        }

        Ok(())
    }
    #[test]
    fn test_order_0_dispatch_boundary() -> io::Result<()> {
        for len in [4095, 4096, 4097] {
            let src = sparse(len, 3, 32);
            for flags in [Flags::empty(), Flags::N32] {
                assert_eq!(decode(&encode(flags, &src)?, len)?, src);
            }
        }
        Ok(())
    }
}
