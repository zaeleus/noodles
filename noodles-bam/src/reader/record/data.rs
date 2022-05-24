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
        fn t(mut src: &[u8], actual: &mut Data, expected: &Data) -> io::Result<()> {
            get_data(&mut src, actual)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = Data::default();

        t(&[], &mut buf, &Data::default())?;

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
