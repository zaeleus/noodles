//! BAM record data field writer.

pub mod field;

use std::io;

use bytes::BufMut;
use noodles_sam::{alignment::record::Data, record::data::field::tag};

use self::field::put_field;

pub(crate) fn put_data<B, D>(dst: &mut B, data: D) -> io::Result<()>
where
    B: BufMut,
    D: Data,
{
    for result in data.iter() {
        let (tag, value) = result?;

        if &tag == tag::CIGAR.as_ref() {
            continue;
        }

        put_field(dst, tag, &value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_sam::{alignment::record_buf::data::field::Value, record::Data as DataBuf};

    use super::*;

    #[test]
    fn test_put_data() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, data: &DataBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_data(buf, data)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let data = DataBuf::default();
        t(&mut buf, &data, &[])?;

        let data = [(tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
            .into_iter()
            .collect();
        t(&mut buf, &data, &[b'N', b'H', b'C', 0x01])?;

        let data = [
            (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
            (tag::READ_GROUP, Value::from("rg0")),
        ]
        .into_iter()
        .collect();
        t(
            &mut buf,
            &data,
            &[
                b'N', b'H', b'C', 0x01, // NH:C:1
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
            ],
        )?;

        Ok(())
    }
}
