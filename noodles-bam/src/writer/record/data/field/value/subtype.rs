use bytes::BufMut;
use noodles_sam::record::data::field::value::Subtype;

pub fn put_subtype<B>(dst: &mut B, subtype: Subtype)
where
    B: BufMut,
{
    dst.put_u8(u8::from(subtype));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_subtype() {
        let mut buf = Vec::new();
        put_subtype(&mut buf, Subtype::Int32);
        assert_eq!(buf, [b'i']);
    }
}
