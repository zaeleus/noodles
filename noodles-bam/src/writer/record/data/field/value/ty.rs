use std::io::{self, Write};

use byteorder::WriteBytesExt;
use noodles_sam::record::data::field::value::Type;

pub fn write_type<W>(writer: &mut W, ty: Type) -> io::Result<()>
where
    W: Write,
{
    writer.write_u8(u8::from(ty))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        let mut buf = Vec::new();
        write_type(&mut buf, Type::Int32)?;
        assert_eq!(buf, [b'i']);
        Ok(())
    }
}
