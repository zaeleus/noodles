use std::{
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::Record;

pub(super) fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    let block_size = u32::try_from(record.block_size())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(block_size)?;

    writer.write_i32::<LittleEndian>(record.ref_id)?;
    writer.write_i32::<LittleEndian>(record.pos)?;

    let l_read_name = u8::try_from(record.read_name.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u8(l_read_name)?;

    let mapq = u8::from(record.mapping_quality());
    writer.write_u8(mapq)?;

    writer.write_u16::<LittleEndian>(record.bin())?;

    let n_cigar_op = u16::try_from(record.cigar().len() / mem::size_of::<u32>())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_cigar_op)?;

    let flag = u16::from(record.flags());
    writer.write_u16::<LittleEndian>(flag)?;

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(l_seq)?;

    writer.write_i32::<LittleEndian>(record.next_ref_id)?;
    writer.write_i32::<LittleEndian>(record.next_pos)?;

    writer.write_i32::<LittleEndian>(record.template_length())?;

    writer.write_all(&record.read_name)?;

    for &raw_op in record.cigar().iter() {
        writer.write_u32::<LittleEndian>(raw_op)?;
    }

    writer.write_all(record.sequence().as_ref())?;
    writer.write_all(record.quality_scores().as_ref())?;
    writer.write_all(record.data().as_ref())?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let record = Record::default();
        write_record(&mut buf, &record)?;

        let expected = [
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

        assert_eq!(buf, expected);

        Ok(())
    }
}
