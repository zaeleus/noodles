mod field;

use std::io;

use self::field::write_field;

pub(super) fn write_typed_map_encoded_data(dst: &mut Vec<u8>, mut src: &[u8]) -> io::Result<()> {
    while !src.is_empty() {
        write_field(dst, &mut src)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_typed_map_encoded_data() -> io::Result<()> {
        let mut dst = Vec::new();

        dst.clear();
        let src = [];
        write_typed_map_encoded_data(&mut dst, &src)?;
        assert_eq!(dst, src);

        dst.clear();
        let src = [b'N', b'H', b'C', 0x01]; // NH:C:1
        write_typed_map_encoded_data(&mut dst, &src)?;
        assert_eq!(dst, src);

        dst.clear();
        let src = [
            b'Z', b'A', b'A', b'n', // ZA:A:n
            b'Z', b'B', b'c', 0x00, // ZB:c:0
            b'Z', b'C', b'C', 0x00, // ZC:C:0
            b'Z', b'D', b's', 0x00, 0x00, // ZD:s:0
            b'Z', b'E', b'S', 0x00, 0x00, // ZE:S:0
            b'Z', b'F', b'i', 0x00, 0x00, 0x00, 0x00, // ZF:i:0
            b'Z', b'G', b'I', 0x00, 0x00, 0x00, 0x00, // ZG:I:0
            b'Z', b'H', b'f', 0x00, 0x00, 0x00, 0x00, // ZH:f:0.0
            b'Z', b'I', b'Z', b'n', b'd', b'l', b's', 0x00, // ZI:Z:ndls
            b'Z', b'J', b'H', b'C', b'A', b'F', b'E', 0x00, // ZJ:H:CAFE
            b'Z', b'K', b'B', b'c', 0x01, 0x00, 0x00, 0x00, 0x00, // ZK:B:c,0
            b'Z', b'L', b'B', b'C', 0x01, 0x00, 0x00, 0x00, 0x00, // ZL:B:C,0
            b'Z', b'M', b'B', b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // ZK:B:s,0
            b'Z', b'N', b'B', b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // ZL:B:S,0
            // ZO:B:i,0
            b'Z', b'O', b'B', b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            // ZP:B:I,0
            b'Z', b'P', b'B', b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            // ZQ:B:f,0
            b'Z', b'Q', b'B', b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        write_typed_map_encoded_data(&mut dst, &src)?;
        assert_eq!(dst, src);

        Ok(())
    }
}
