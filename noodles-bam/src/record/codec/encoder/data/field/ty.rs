use noodles_sam::alignment::record::data::field::Type;

use crate::record::codec::encoder::num::write_u8;

pub(super) fn write_type(dst: &mut Vec<u8>, ty: Type) {
    let n = encode(ty);
    write_u8(dst, n);
}

pub fn encode(ty: Type) -> u8 {
    match ty {
        Type::Character => b'A',
        Type::Int8 => b'c',
        Type::UInt8 => b'C',
        Type::Int16 => b's',
        Type::UInt16 => b'S',
        Type::Int32 => b'i',
        Type::UInt32 => b'I',
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
    fn test_write_type() {
        let mut buf = Vec::new();
        write_type(&mut buf, Type::Character);
        assert_eq!(buf, [b'A']);
    }

    #[test]
    fn test_type_to_u8() {
        assert_eq!(encode(Type::Character), b'A');
        assert_eq!(encode(Type::Int8), b'c');
        assert_eq!(encode(Type::UInt8), b'C');
        assert_eq!(encode(Type::Int16), b's');
        assert_eq!(encode(Type::UInt16), b'S');
        assert_eq!(encode(Type::Int32), b'i');
        assert_eq!(encode(Type::UInt32), b'I');
        assert_eq!(encode(Type::Float), b'f');
        assert_eq!(encode(Type::String), b'Z');
        assert_eq!(encode(Type::Hex), b'H');
        assert_eq!(encode(Type::Array), b'B');
    }
}
