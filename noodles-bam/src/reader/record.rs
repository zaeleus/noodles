use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::{self as sam, alignment::Record};

pub(crate) fn read_record<R>(
    reader: &mut R,
    header: &sam::Header,
    buf: &mut Vec<u8>,
    record: &mut Record,
) -> io::Result<usize>
where
    R: Read,
{
    use crate::record::codec::decode;

    let block_size = match read_raw_record(reader, buf)? {
        0 => return Ok(0),
        n => n,
    };

    let mut src = &buf[..];
    decode(&mut src, header, record).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(block_size)
}

pub(super) fn read_raw_record<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: Read,
{
    let block_size = match read_block_size(reader)? {
        0 => return Ok(0),
        n => n,
    };

    buf.resize(block_size, 0);
    reader.read_exact(buf)?;

    Ok(block_size)
}

fn read_block_size<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read,
{
    match reader.read_u32::<LittleEndian>() {
        Ok(n) => usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(0),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_block_size() -> io::Result<()> {
        let data = [0x08, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader)?, 8);

        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader)?, 0);

        Ok(())
    }

    #[test]
    fn test_read_record() -> io::Result<()> {
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
        let header = sam::Header::default();
        let mut buf = Vec::new();
        let mut record = Record::default();
        let block_size = read_record(&mut reader, &header, &mut buf, &mut record)?;

        assert_eq!(block_size, 34);
        assert_eq!(record, Record::default());

        Ok(())
    }
}
