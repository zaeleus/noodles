//! BAM record data field component writers.

mod tag;
mod value;

use std::io;

use bytes::BufMut;
use noodles_sam::record::data::Field;

use self::tag::put_tag;
pub use self::value::put_value;

pub(super) fn put_field<B>(dst: &mut B, field: &Field) -> io::Result<()>
where
    B: BufMut,
{
    put_tag(dst, field.tag());
    value::put_type(dst, field.value().ty());
    put_value(dst, field.value())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_field() -> io::Result<()> {
        use noodles_sam::record::data::field::{Tag, Value};

        let mut buf = Vec::new();
        let field = Field::new(Tag::AlignmentHitCount, Value::from(1));
        put_field(&mut buf, &field)?;
        assert_eq!(buf, [b'N', b'H', b'C', 0x01]);

        Ok(())
    }
}
