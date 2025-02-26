//! BAM record data field writer.

pub mod field;

use std::io;

use noodles_sam::alignment::record::{data::field::Tag, Data};

use self::field::write_field;

pub(super) fn write_data<D>(dst: &mut Vec<u8>, data: D) -> io::Result<()>
where
    D: Data,
{
    for result in data.iter() {
        let (tag, value) = result?;

        if &tag == Tag::CIGAR.as_ref() {
            continue;
        }

        write_field(dst, tag, &value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record_buf::{data::field::Value, Data as DataBuf};

    use super::*;

    #[test]
    fn test_write_data() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, data: &DataBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data(buf, data)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let data = DataBuf::default();
        t(&mut buf, &data, &[])?;

        let data = [(Tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
            .into_iter()
            .collect();
        t(&mut buf, &data, &[b'N', b'H', b'C', 0x01])?;

        let data = [
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
            (Tag::READ_GROUP, Value::from("rg0")),
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
