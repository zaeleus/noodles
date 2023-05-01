//! BAM record data field writer.

pub mod field;

use std::io;

use bytes::BufMut;
use noodles_sam::record::{data::field::tag, Data};

use self::field::put_field;

pub(crate) fn put_data<B>(dst: &mut B, data: &Data) -> io::Result<()>
where
    B: BufMut,
{
    for (tag, value) in data.iter() {
        if tag == tag::CIGAR {
            continue;
        }

        put_field(dst, tag, value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_data() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, data: &Data, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_data(buf, data)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Data::default(), &[])?;
        t(&mut buf, &"NH:i:1".parse()?, &[b'N', b'H', b'C', 0x01])?;
        t(
            &mut buf,
            &"NH:i:1\tRG:Z:rg0".parse()?,
            &[
                b'N', b'H', b'C', 0x01, // NH:C:1
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
            ],
        )?;

        Ok(())
    }
}
