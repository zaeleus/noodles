use std::io;

use crate::alignment::record::data::field::Value;

pub(super) fn parse_int32_value<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    let (n, i) = lexical_core::parse_partial(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *src = &src[i..];

    Ok(Value::Int32(n))
}
