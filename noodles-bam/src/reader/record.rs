//! BAM record field readers.

pub mod data;

use std::{
    io::{self, Read},
    mem,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam as sam;

use crate::{
    record::{Cigar, Data, QualityScores, Sequence},
    Record,
};

pub(crate) fn read_record<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    record: &mut Record,
) -> io::Result<usize>
where
    R: Read,
{
    let block_size = match reader.read_u32::<LittleEndian>() {
        Ok(bs) => usize::try_from(bs).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    };

    buf.resize(block_size, Default::default());
    reader.read_exact(buf)?;

    let mut reader = &buf[..];
    let reader = &mut reader;

    record.ref_id = reader.read_i32::<LittleEndian>()?;
    record.pos = reader.read_i32::<LittleEndian>()?;

    let l_read_name = reader.read_u8().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    *record.mapping_quality_mut() = read_mapping_quality(reader)?;
    *record.bin_mut() = reader.read_u16::<LittleEndian>()?;

    let n_cigar_op = reader.read_u16::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    *record.flags_mut() = read_flag(reader)?;

    let l_seq = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    record.next_ref_id = reader.read_i32::<LittleEndian>()?;
    record.next_pos = reader.read_i32::<LittleEndian>()?;

    *record.template_length_mut() = reader.read_i32::<LittleEndian>()?;

    read_read_name(reader, &mut record.read_name, l_read_name)?;
    read_cigar(reader, record.cigar_mut(), n_cigar_op)?;
    read_seq(reader, record.sequence_mut(), l_seq)?;
    read_qual(reader, record.quality_scores_mut(), l_seq)?;
    read_data(
        reader,
        record.data_mut(),
        block_size,
        l_read_name,
        n_cigar_op,
        l_seq,
    )?;

    Ok(block_size)
}

fn read_mapping_quality<R>(reader: &mut R) -> io::Result<sam::record::MappingQuality>
where
    R: Read,
{
    reader.read_u8().map(sam::record::MappingQuality::from)
}

fn read_flag<R>(reader: &mut R) -> io::Result<sam::record::Flags>
where
    R: Read,
{
    reader
        .read_u16::<LittleEndian>()
        .map(sam::record::Flags::from)
}

fn read_read_name<R>(reader: &mut R, read_name: &mut Vec<u8>, l_read_name: usize) -> io::Result<()>
where
    R: Read,
{
    read_name.resize(l_read_name, Default::default());
    reader.read_exact(read_name)?;
    Ok(())
}

fn read_cigar<R>(reader: &mut R, cigar: &mut Cigar, n_cigar_op: usize) -> io::Result<()>
where
    R: Read,
{
    cigar.resize(n_cigar_op, Default::default());
    reader.read_u32_into::<LittleEndian>(cigar)?;
    Ok(())
}

fn read_seq<R>(reader: &mut R, sequence: &mut Sequence, l_seq: usize) -> io::Result<()>
where
    R: Read,
{
    sequence.set_len(l_seq);

    let seq = sequence.as_mut();
    let seq_len = (l_seq + 1) / 2;
    seq.resize(seq_len, Default::default());
    reader.read_exact(seq)?;

    Ok(())
}

fn read_qual<R>(reader: &mut R, qual: &mut QualityScores, l_seq: usize) -> io::Result<()>
where
    R: Read,
{
    qual.resize(l_seq, Default::default());
    reader.read_exact(qual)?;
    Ok(())
}

fn read_data<R>(
    reader: &mut R,
    data: &mut Data,
    block_size: usize,
    l_read_name: usize,
    n_cigar_op: usize,
    l_seq: usize,
) -> io::Result<()>
where
    R: Read,
{
    let cigar_len = mem::size_of::<u32>() * n_cigar_op;
    let seq_len = (l_seq + 1) / 2;
    let data_offset = 32 + l_read_name + cigar_len + seq_len + l_seq;
    let data_len = block_size - data_offset;

    let buf = data.as_mut();
    buf.resize(data_len, Default::default());
    reader.read_exact(buf)?;

    data.index()?;

    Ok(())
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
        let block_size = read_record(&mut reader, &mut Vec::new(), &mut record)?;

        assert_eq!(block_size, 34);
        assert_eq!(record, Record::default());

        Ok(())
    }
}
