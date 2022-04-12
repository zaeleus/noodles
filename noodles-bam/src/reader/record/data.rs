//! BAM record data field reader.

pub mod field;

pub(crate) use self::field::get_field;

use std::io;

use bytes::Buf;
use noodles_sam as sam;

pub(super) fn get_data<B>(src: &mut B, data: &mut sam::record::Data) -> io::Result<()>
where
    B: Buf,
{
    data.clear();

    while let Some(field) = get_field(src)? {
        data.insert(field);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_data() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            mut src: &[u8],
            actual: &mut sam::record::Data,
            expected: &sam::record::Data,
        ) -> io::Result<()> {
            get_data(&mut src, actual)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = sam::record::Data::default();

        t(&[], &mut buf, &sam::record::Data::default())?;

        t(
            &[b'N', b'H', b'C', 0x01], // NH:C:0
            &mut buf,
            &"NH:i:1".parse()?,
        )?;

        t(
            &[
                b'N', b'H', b'C', 0x01, // NH:C:0
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
            ],
            &mut buf,
            &"NH:i:1\tRG:Z:rg0".parse()?,
        )?;

        Ok(())
    }
}
