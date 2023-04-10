//! BAM record data field reader.

pub mod field;

pub(crate) use self::field::get_field;

use std::io;

use bytes::Buf;
use noodles_sam::record::Data;

pub(crate) fn get_data<B>(src: &mut B, data: &mut Data) -> io::Result<()>
where
    B: Buf,
{
    data.clear();

    while let Some((tag, value)) = get_field(src)? {
        data.insert(tag, value);
    }

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_get_data() -> io::Result<()> {
        use noodles_sam::record::data::field::{Tag, Value};

        fn t(mut src: &[u8], actual: &mut Data, expected: &Data) -> io::Result<()> {
            get_data(&mut src, actual)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = Data::default();

        let expected = Data::default();
        t(&[], &mut buf, &expected)?;

        let expected = [(Tag::AlignmentHitCount, Value::UInt8(1))]
            .into_iter()
            .collect();

        t(
            &[b'N', b'H', b'C', 0x01], // NH:C:0
            &mut buf,
            &expected,
        )?;

        let expected = [
            (Tag::AlignmentHitCount, Value::UInt8(1)),
            (Tag::ReadGroup, Value::String(String::from("rg0"))),
        ]
        .into_iter()
        .collect();

        t(
            &[
                b'N', b'H', b'C', 0x01, // NH:C:0
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
            ],
            &mut buf,
            &expected,
        )?;

        Ok(())
    }
}
