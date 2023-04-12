//! BAM record data field component readers.

mod tag;
mod value;

pub use self::value::get_value;

use std::io;

use bytes::Buf;
use noodles_sam::record::data::field::{Tag, Value};

pub(crate) fn get_field<B>(src: &mut B) -> io::Result<(Tag, Value)>
where
    B: Buf,
{
    use self::tag::get_tag;

    let tag = get_tag(src)?;

    let ty = value::get_type(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let value = get_value(src, ty).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((tag, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_field() -> io::Result<()> {
        use noodles_sam::record::data::field::{Tag, Value};

        let data = [b'N', b'H', b'C', 0x01];
        let mut reader = &data[..];
        let actual = get_field(&mut reader)?;
        let expected = (Tag::AlignmentHitCount, Value::from(1));
        assert_eq!(actual, expected);

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
           get_field(&mut reader),
           Err(e) if e.kind() == io::ErrorKind::UnexpectedEof,
        ));

        Ok(())
    }
}
