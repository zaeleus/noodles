use std::io::{self, Write};

use crate::record::data::field::value::Type;

pub fn write_type<W>(writer: &mut W, mut ty: Type) -> io::Result<()>
where
    W: Write,
{
    if matches!(
        ty,
        Type::Int8 | Type::UInt8 | Type::Int16 | Type::UInt16 | Type::Int32 | Type::UInt32
    ) {
        ty = Type::Int32;
    }

    let c = u8::from(ty);
    writer.write_all(&[c])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, ty: Type, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_type(buf, ty)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Type::Character, b"A")?;
        t(&mut buf, Type::Int8, b"i")?;
        t(&mut buf, Type::UInt8, b"i")?;
        t(&mut buf, Type::Int16, b"i")?;
        t(&mut buf, Type::UInt16, b"i")?;
        t(&mut buf, Type::Int32, b"i")?;
        t(&mut buf, Type::UInt32, b"i")?;
        t(&mut buf, Type::Float, b"f")?;
        t(&mut buf, Type::String, b"Z")?;
        t(&mut buf, Type::Hex, b"H")?;
        t(&mut buf, Type::Array, b"B")?;

        Ok(())
    }
}
