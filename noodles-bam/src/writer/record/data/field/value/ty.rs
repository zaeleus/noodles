use bytes::BufMut;
use noodles_sam::record::data::field::value::Type;

pub fn put_type<B>(dst: &mut B, ty: Type)
where
    B: BufMut,
{
    dst.put_u8(u8::from(ty));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_type() {
        let mut buf = Vec::new();
        put_type(&mut buf, Type::Int32);
        assert_eq!(buf, [b'i']);
    }
}
