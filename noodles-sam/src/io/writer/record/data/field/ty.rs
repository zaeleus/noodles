use std::io::{self, Write};

use crate::alignment::record::data::field::Type;

pub fn write_type<W>(writer: &mut W, ty: Type) -> io::Result<()>
where
    W: Write,
{
    let c = encode(ty);
    writer.write_all(&[c])
}

// § 1.5 "The alignment section: optional fields" (2024-11-06): "In an optional field, `TYPE` is a
// single case-sensitive letter which defines the format of `VALUE`: `[AifZHB]`."
fn encode(ty: Type) -> u8 {
    match ty {
        Type::Character => b'A',
        Type::Int8 | Type::UInt8 | Type::Int16 | Type::UInt16 | Type::Int32 | Type::UInt32 => b'i',
        Type::Float => b'f',
        Type::String => b'Z',
        Type::Hex => b'H',
        Type::Array => b'B',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        let mut buf = Vec::new();
        write_type(&mut buf, Type::Character)?;
        assert_eq!(buf, b"A");
        Ok(())
    }

    #[test]
    fn test_encode() {
        assert_eq!(encode(Type::Character), b'A');
        assert_eq!(encode(Type::Int8), b'i');
        assert_eq!(encode(Type::UInt8), b'i');
        assert_eq!(encode(Type::Int16), b'i');
        assert_eq!(encode(Type::UInt16), b'i');
        assert_eq!(encode(Type::Int32), b'i');
        assert_eq!(encode(Type::UInt32), b'i');
        assert_eq!(encode(Type::Float), b'f');
        assert_eq!(encode(Type::String), b'Z');
        assert_eq!(encode(Type::Hex), b'H');
        assert_eq!(encode(Type::Array), b'B');
    }
}
