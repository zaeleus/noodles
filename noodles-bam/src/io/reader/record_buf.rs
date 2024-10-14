use std::io::{self, Read};

use noodles_sam::alignment::RecordBuf;

use super::read_record;

pub(crate) fn read_record_buf<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    record: &mut RecordBuf,
) -> io::Result<usize>
where
    R: Read,
{
    use crate::record::codec::decode;

    let block_size = match read_record(reader, buf)? {
        0 => return Ok(0),
        n => n,
    };

    let mut src = &buf[..];
    decode(&mut src, record).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(block_size)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_record_buf() -> io::Result<()> {
        let data = [
            0x22, 0x00, 0x00, 0x00, // block_size = 34
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x00, 0x00, // n_cigar_op = 0
            0x04, 0x00, // flag = 4
            0x00, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            0x2a, 0x00, // read_name = "*\x00"
        ];

        let mut reader = &data[..];
        let mut buf = Vec::new();
        let mut record = RecordBuf::default();
        let block_size = read_record_buf(&mut reader, &mut buf, &mut record)?;

        assert_eq!(block_size, 34);
        assert_eq!(record, RecordBuf::default());

        Ok(())
    }
}
