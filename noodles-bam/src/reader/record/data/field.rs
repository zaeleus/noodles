//! BAM record data field component readers.

mod tag;
mod value;

pub use self::value::get_value;

use std::io;

use bytes::Buf;
use noodles_sam::record::data::field::{Tag, Value};

pub(crate) fn get_field<B>(src: &mut B) -> io::Result<Option<(Tag, Value)>>
where
    B: Buf,
{
    use self::tag::get_tag;

    let tag = match get_tag(src) {
        Ok(t) => t,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    };

    let ty = value::get_type(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let value = get_value(src, ty).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(Some((tag, value)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_field() -> io::Result<()> {
        use noodles_sam::record::data::field::{Tag, Value};

        let data = [];
        let mut reader = &data[..];
        assert!(get_field(&mut reader)?.is_none());

        let data = [b'N', b'H', b'C', 0x01];
        let mut reader = &data[..];
        let actual = get_field(&mut reader)?;
        let expected = (Tag::AlignmentHitCount, Value::from(1));
        assert_eq!(actual, Some(expected));

        Ok(())
    }
}
