//! BAM record data field component writers.

mod tag;
mod value;

use std::io;

use bytes::BufMut;
use noodles_sam::record::data::field::{Tag, Value};

use self::tag::put_tag;
pub use self::value::put_value;

pub(super) fn put_field<B>(dst: &mut B, tag: Tag, value: &Value) -> io::Result<()>
where
    B: BufMut,
{
    put_tag(dst, tag);
    value::put_type(dst, value.ty());
    put_value(dst, value)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_field() -> io::Result<()> {
        let mut buf = Vec::new();
        let (tag, value) = (Tag::AlignmentHitCount, Value::from(1));
        put_field(&mut buf, tag, &value)?;
        assert_eq!(buf, [b'N', b'H', b'C', 0x01]);
        Ok(())
    }
}
