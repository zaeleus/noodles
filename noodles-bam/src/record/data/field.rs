//! BAM record data field.

mod tag;
mod ty;
pub mod value;

use std::io;

use noodles_sam::alignment::record::data::field::{Tag, Value};

pub(crate) use self::{tag::decode_tag, ty::decode_type, value::decode_value};

pub(super) fn decode_field<'a>(src: &mut &'a [u8]) -> io::Result<(Tag, Value<'a>)> {
    let tag = decode_tag(src)?;

    let ty = decode_type(src)?;
    let value = decode_value(src, ty)?;

    Ok((tag, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_field() -> io::Result<()> {
        let mut src = &[b'N', b'H', b'C', 0x01][..];
        assert!(matches!(
            decode_field(&mut src)?,
            (Tag::ALIGNMENT_HIT_COUNT, Value::UInt8(1))
        ));

        let mut src = &[][..];
        assert!(matches!(
            decode_field(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }
}
