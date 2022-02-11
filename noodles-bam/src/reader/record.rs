//! BAM record field readers.

pub mod data;

use std::{
    io::{self, Read},
    mem,
    num::NonZeroUsize,
};

use byteorder::{LittleEndian, ReadBytesExt};
use bytes::Buf;
use noodles_sam as sam;

use crate::{
    record::{Cigar, Data, QualityScores, ReferenceSequenceId, Sequence},
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

    read_record_buf(&buf[..], record)?;

    Ok(block_size)
}

pub(crate) fn read_record_buf<B>(mut buf: B, record: &mut Record) -> io::Result<()>
where
    B: Buf,
{
    *record.reference_sequence_id_mut() = read_reference_sequence_id(&mut buf)?;
    record.pos = read_position(&mut buf)?;

    let l_read_name = NonZeroUsize::new(usize::from(buf.get_u8()))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid l_read_name"))?;

    *record.mapping_quality_mut() = read_mapping_quality(&mut buf)?;
    *record.bin_mut() = buf.get_u16_le();

    let n_cigar_op = usize::from(buf.get_u16_le());

    *record.flags_mut() = read_flag(&mut buf)?;

    let l_seq = usize::try_from(buf.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *record.mate_reference_sequence_id_mut() = read_reference_sequence_id(&mut buf)?;
    record.next_pos = read_position(&mut buf)?;

    *record.template_length_mut() = buf.get_i32_le();

    read_read_name(&mut buf, &mut record.read_name, l_read_name)?;
    read_cigar(&mut buf, record.cigar_mut(), n_cigar_op)?;
    read_seq(&mut buf, record.sequence_mut(), l_seq)?;
    read_qual(&mut buf, record.quality_scores_mut(), l_seq)?;

    read_data(&mut buf, record.data_mut())?;

    Ok(())
}

fn read_reference_sequence_id<B>(buf: &mut B) -> io::Result<Option<ReferenceSequenceId>>
where
    B: Buf,
{
    use crate::record::reference_sequence_id::UNMAPPED;

    if buf.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    match buf.get_i32_le() {
        UNMAPPED => Ok(None),
        n => ReferenceSequenceId::try_from(n)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

fn read_position<B>(buf: &mut B) -> io::Result<i32>
where
    B: Buf,
{
    if buf.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    Ok(buf.get_i32_le())
}

fn read_mapping_quality<B>(buf: &mut B) -> io::Result<Option<sam::record::MappingQuality>>
where
    B: Buf,
{
    use sam::record::mapping_quality::MISSING;

    if buf.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    match buf.get_u8() {
        MISSING => Ok(None),
        n => sam::record::MappingQuality::try_from(n)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

fn read_flag<B>(buf: &mut B) -> io::Result<sam::record::Flags>
where
    B: Buf,
{
    if buf.remaining() < mem::size_of::<u16>() {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    Ok(sam::record::Flags::from(buf.get_u16_le()))
}

fn read_read_name<B>(
    buf: &mut B,
    read_name: &mut Vec<u8>,
    l_read_name: NonZeroUsize,
) -> io::Result<()>
where
    B: Buf,
{
    let len = usize::from(l_read_name);

    if buf.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    read_name.resize(len, Default::default());
    buf.copy_to_slice(read_name);

    Ok(())
}

fn read_cigar<B>(buf: &mut B, cigar: &mut Cigar, n_cigar_op: usize) -> io::Result<()>
where
    B: Buf,
{
    if buf.remaining() < mem::size_of::<u32>() * n_cigar_op {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    let cigar = cigar.as_mut();

    cigar.resize(n_cigar_op, Default::default());
    cigar.clear();

    for _ in 0..n_cigar_op {
        cigar.push(buf.get_u32_le());
    }

    Ok(())
}

fn read_seq<B>(buf: &mut B, sequence: &mut Sequence, l_seq: usize) -> io::Result<()>
where
    B: Buf,
{
    sequence.set_len(l_seq);

    let seq = sequence.as_mut();
    let seq_len = (l_seq + 1) / 2;

    if buf.remaining() < seq_len {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    seq.resize(seq_len, Default::default());
    buf.copy_to_slice(seq);

    Ok(())
}

fn read_qual<B>(buf: &mut B, quality_scores: &mut QualityScores, l_seq: usize) -> io::Result<()>
where
    B: Buf,
{
    if buf.remaining() < l_seq {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    let qual = quality_scores.as_mut();
    qual.resize(l_seq, Default::default());
    buf.copy_to_slice(qual);

    Ok(())
}

fn read_data<B>(buf: &mut B, data: &mut Data) -> io::Result<()>
where
    B: Buf,
{
    let data_len = buf.remaining();

    let data_buf = data.as_mut();
    data_buf.resize(data_len, Default::default());
    buf.copy_to_slice(data_buf);

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

    #[test]
    fn test_read_record_buf_with_invalid_l_read_name() {
        let data = [
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x00, // l_read_name = 0
        ];

        let mut reader = &data[..];
        let mut record = Record::default();

        assert!(matches!(
            read_record_buf(&mut reader, &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
