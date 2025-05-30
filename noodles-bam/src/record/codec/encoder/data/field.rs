//! BAM record data field component writers.

mod tag;
pub mod ty;
mod value;

use std::io;

use noodles_sam::alignment::record::{
    Cigar,
    data::field::{Tag, Type, Value, value::array::Subtype},
};

pub use self::value::write_value;
use self::{tag::write_tag, ty::write_type};

pub(super) fn write_field(dst: &mut Vec<u8>, tag: Tag, value: &Value) -> io::Result<()> {
    write_tag(dst, tag);
    write_type(dst, value.ty());
    write_value(dst, value)?;
    Ok(())
}

pub(crate) fn write_cigar<C>(dst: &mut Vec<u8>, cigar: &C) -> io::Result<()>
where
    C: Cigar,
{
    write_tag(dst, Tag::CIGAR);
    write_type(dst, Type::Array);
    value::array::write_header(dst, Subtype::UInt32, cigar.len())?;
    crate::record::codec::encoder::write_cigar(dst, cigar)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let (tag, value) = (Tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1));
        write_field(&mut buf, tag, &value)?;
        assert_eq!(buf, [b'N', b'H', b'C', 0x01]);
        Ok(())
    }

    #[test]
    fn test_write_cigar() -> io::Result<()> {
        use noodles_sam::alignment::{
            record::cigar::{Op, op::Kind},
            record_buf::Cigar as CigarBuf,
        };

        let mut buf = Vec::new();
        let cigar: CigarBuf = [Op::new(Kind::Match, 4)].into_iter().collect();
        write_cigar(&mut buf, &cigar)?;

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
