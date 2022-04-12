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
