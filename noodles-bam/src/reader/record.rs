//! BAM record field readers.

pub mod data;

use std::{
    io::{self, Read},
    mem,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam as sam;

use crate::Record;

pub(crate) fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: Read,
{
    let block_size = match reader.read_u32::<LittleEndian>() {
        Ok(bs) => usize::try_from(bs).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    };

    record.ref_id = reader.read_i32::<LittleEndian>()?;
    record.pos = reader.read_i32::<LittleEndian>()?;

    let l_read_name = reader.read_u8().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    record.mapq = reader.read_u8()?;
    record.bin = reader.read_u16::<LittleEndian>()?;

    let n_cigar_op = reader.read_u16::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    *record.flags_mut() = read_flag(reader)?;

    let l_seq = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    record.next_ref_id = reader.read_i32::<LittleEndian>()?;
    record.next_pos = reader.read_i32::<LittleEndian>()?;
    record.tlen = reader.read_i32::<LittleEndian>()?;

    record.read_name.resize(l_read_name, Default::default());
    reader.read_exact(&mut record.read_name)?;

    let cigar = record.cigar_mut();
    cigar.resize(n_cigar_op, Default::default());
    reader.read_u32_into::<LittleEndian>(cigar)?;

    let seq_len = (l_seq + 1) / 2;

    let seq = record.sequence_mut();
    seq.set_base_count(l_seq);

    seq.resize(seq_len, Default::default());
    reader.read_exact(seq)?;

    let qual = record.quality_scores_mut();
    qual.resize(l_seq, Default::default());
    reader.read_exact(qual)?;

    let data = record.data_mut();

    let cigar_len = mem::size_of::<u32>() * n_cigar_op;
    let data_offset = 32 + l_read_name + cigar_len + seq_len + l_seq;
    let data_len = block_size - data_offset;
    data.resize(data_len, Default::default());

    reader.read_exact(data)?;

    Ok(block_size)
}

fn read_flag<R>(reader: &mut R) -> io::Result<sam::record::Flags>
where
    R: Read,
{
    reader
        .read_u16::<LittleEndian>()
        .map(sam::record::Flags::from)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let mut record = Record::default();
        let block_size = read_record(&mut reader, &mut record)?;

        assert_eq!(block_size, 34);
        assert_eq!(record, Record::default());

        Ok(())
    }
}
