//! BAM record data field component writers.

mod tag;
pub mod ty;
mod value;

use std::io;

use bytes::BufMut;
use noodles_sam::alignment::record::{
    data::field::{value::array::Subtype, Tag, Type, Value},
    fields::Cigar,
};

pub use self::value::put_value;
use self::{tag::put_tag, ty::put_type};

pub(super) fn put_field<B>(dst: &mut B, tag: Tag, value: &Value) -> io::Result<()>
where
    B: BufMut,
{
    put_tag(dst, tag);
    put_type(dst, value.ty());
    put_value(dst, value)?;
    Ok(())
}

pub(crate) fn put_cigar<B, C>(dst: &mut B, cigar: &C) -> io::Result<()>
where
    B: BufMut,
    C: Cigar,
{
    put_tag(dst, Tag::CIGAR);
    put_type(dst, Type::Array);
    value::array::put_header(dst, Subtype::UInt32, cigar.len())?;
    crate::record::codec::encoder::put_cigar(dst, cigar)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let (tag, value) = (Tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1));
        put_field(&mut buf, tag, &value)?;
        assert_eq!(buf, [b'N', b'H', b'C', 0x01]);
        Ok(())
    }

    #[test]
    fn test_put_cigar() -> io::Result<()> {
        use noodles_sam::alignment::{
            record::cigar::{op::Kind, Op},
            record_buf::Cigar as CigarBuf,
        };

        let mut buf = Vec::new();
        let cigar: CigarBuf = [Op::new(Kind::Match, 4)].into_iter().collect();
        put_cigar(&mut buf, &cigar)?;

        let expected = [
            b'C', b'G', // tag = CG
            b'B', b'I', // type = [u32]
            0x01, 0x00, 0x00, 0x00, // count = 1
            0x40, 0x00, 0x00, 0x00, // cigar = 4M
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
