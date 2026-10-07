mod field;

use std::io::{self, Write};

use self::field::write_field;

pub(super) fn write_typed_map_encoded_data<W>(writer: &mut W, mut src: &[u8]) -> io::Result<()>
where
    W: Write,
{
    while !src.is_empty() {
        write_field(writer, &mut src)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_typed_map_encoded_data() -> io::Result<()> {
        let mut buf = Vec::new();

        let src = [
            b'N', b'H', b'c', 0x01, // NH:c:1
            b'C', b'O', b'Z', b'n', b'd', b'l', b's', 0x00, // CO:Z:ndls
        ];

        write_typed_map_encoded_data(&mut buf, &src)?;
        assert_eq!(buf, b"\tNH:i:1\tCO:Z:ndls");

        Ok(())
    }
}
