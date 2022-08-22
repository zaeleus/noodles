use std::io::{self, Write};

use byteorder::WriteBytesExt;

use super::Flags;
use crate::writer::num::write_uint7;

#[allow(dead_code)]
pub fn rans_encode_nx16(flags: Flags, src: &[u8]) -> io::Result<Vec<u8>> {
    let mut dst = Vec::new();

    dst.write_u8(u8::from(flags))?;

    if !flags.contains(Flags::NO_SIZE) {
        let n =
            u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7(&mut dst, n)?;
    }

    let _n = if flags.contains(Flags::N32) { 32 } else { 4 };

    if flags.contains(Flags::STRIPE) {
        todo!("rans_encode_stripe");
    }

    if flags.contains(Flags::PACK) {
        todo!("encode_pack_meta");
    }

    if flags.contains(Flags::RLE) {
        todo!("encode_rle_meta");
    }

    if flags.contains(Flags::CAT) {
        dst.write_all(src)?;
    } else if flags.contains(Flags::ORDER) {
        todo!("rans_encode_nx16_1");
    } else {
        todo!("rans_encode_nx16_0");
    }

    Ok(dst)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rans_encode_nx16_uncompressed() -> io::Result<()> {
        let actual = rans_encode_nx16(Flags::CAT, b"noodles")?;

        let expected = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
