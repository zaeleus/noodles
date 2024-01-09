use bytes::BufMut;
use noodles_sam::record::data::field::Type;

pub(super) fn put_type<B>(dst: &mut B, ty: Type)
where
    B: BufMut,
{
    let n = type_to_u8(ty);
    dst.put_u8(n);
}

pub fn type_to_u8(ty: Type) -> u8 {
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
    fn test_put_type() {
        let mut buf = Vec::new();
        put_type(&mut buf, Type::Character);
        assert_eq!(buf, [b'A']);
    }

    #[test]
    fn test_type_to_u8() {
        assert_eq!(type_to_u8(Type::Character), b'A');
        assert_eq!(type_to_u8(Type::Int8), b'c');
        assert_eq!(type_to_u8(Type::UInt8), b'C');
        assert_eq!(type_to_u8(Type::Int16), b's');
        assert_eq!(type_to_u8(Type::UInt16), b'S');
        assert_eq!(type_to_u8(Type::Int32), b'i');
        assert_eq!(type_to_u8(Type::UInt32), b'I');
        assert_eq!(type_to_u8(Type::Float), b'f');
        assert_eq!(type_to_u8(Type::String), b'Z');
        assert_eq!(type_to_u8(Type::Hex), b'H');
        assert_eq!(type_to_u8(Type::Array), b'B');
    }
}
