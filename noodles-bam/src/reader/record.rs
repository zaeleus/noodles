//! BAM record field readers.

mod cigar;
pub mod data;
mod mapping_quality;
mod quality_scores;
mod read_name;
mod sequence;

use std::{
    io::{self, Read},
    mem,
    num::NonZeroUsize,
};

use byteorder::{LittleEndian, ReadBytesExt};
use bytes::Buf;
use noodles_core::Position;
use noodles_sam::{self as sam, reader::record::Fields};

use crate::Record;

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

    let mut src = &buf[..];
    decode_record(&mut src, record)?;

    Ok(block_size)
}

pub(super) fn read_record_with_fields<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    record: &mut Record,
    fields: Fields,
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

    let mut src = &buf[..];
    decode_record_with_fields(&mut src, record, fields)?;

    Ok(block_size)
}

pub(crate) fn decode_record<B>(src: &mut B, record: &mut Record) -> io::Result<()>
where
    B: Buf,
{
    decode_record_with_fields(src, record, Fields::all())
}

pub(crate) fn decode_record_with_fields<B>(
    src: &mut B,
    record: &mut Record,
    fields: Fields,
) -> io::Result<()>
where
    B: Buf,
{
    use self::{
        cigar::get_cigar, data::get_data, mapping_quality::get_mapping_quality,
        quality_scores::get_quality_scores, read_name::get_read_name, sequence::get_sequence,
    };

    if fields.contains(Fields::REFERENCE_SEQUENCE_ID) {
        *record.reference_sequence_id_mut() = get_reference_sequence_id(src)?;
    } else {
        src.advance(mem::size_of::<i32>());
    }

    if fields.contains(Fields::ALIGNMENT_START) {
        *record.position_mut() = get_position(src)?;
    } else {
        src.advance(mem::size_of::<i32>());
    }

    let l_read_name = NonZeroUsize::new(usize::from(src.get_u8()))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid l_read_name"))?;

    if fields.contains(Fields::MAPPING_QUALITY) {
        *record.mapping_quality_mut() = get_mapping_quality(src)?;
    } else {
        src.advance(mem::size_of::<u8>());
    }

    // Discard bin.
    src.advance(mem::size_of::<u16>());

    let n_cigar_op = usize::from(src.get_u16_le());

    if fields.contains(Fields::FLAGS) {
        *record.flags_mut() = get_flags(src)?;
    } else {
        src.advance(mem::size_of::<u16>());
    }

    let l_seq = usize::try_from(src.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    if fields.contains(Fields::MATE_REFERENCE_SEQUENCE_ID) {
        *record.mate_reference_sequence_id_mut() = get_reference_sequence_id(src)?;
    } else {
        src.advance(mem::size_of::<i32>());
    }

    if fields.contains(Fields::MATE_ALIGNMENT_START) {
        *record.mate_position_mut() = get_position(src)?;
    } else {
        src.advance(mem::size_of::<i32>());
    }

    if fields.contains(Fields::TEMPLATE_LENGTH) {
        *record.template_length_mut() = src.get_i32_le();
    } else {
        src.advance(mem::size_of::<i32>());
    }

    if fields.contains(Fields::READ_NAME) {
        get_read_name(src, record.read_name_mut(), l_read_name)?;
    } else {
        src.advance(usize::from(l_read_name));
    }

    if fields.contains(Fields::CIGAR) {
        get_cigar(src, record.cigar_mut(), n_cigar_op)?;
    } else {
        let len = mem::size_of::<u32>() * n_cigar_op;
        src.advance(len);
    }

    if fields.contains(Fields::SEQUENCE) {
        get_sequence(src, record.sequence_mut(), l_seq)?;
    } else {
        let len = (l_seq + 1) / 2;
        src.advance(len);
    }

    if fields.contains(Fields::QUALITY_SCORES) {
        get_quality_scores(src, record.quality_scores_mut(), l_seq)?;
    } else {
        src.advance(l_seq);
    }

    if fields.contains(Fields::DATA) {
        get_data(src, record.data_mut())?;
    }

    Ok(())
}

fn get_reference_sequence_id<B>(src: &mut B) -> io::Result<Option<usize>>
where
    B: Buf,
{
    const UNMAPPED: i32 = -1;

    if src.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_i32_le() {
        UNMAPPED => Ok(None),
        n => usize::try_from(n)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

fn get_position<B>(src: &mut B) -> io::Result<Option<Position>>
where
    B: Buf,
{
    use crate::record::UNMAPPED_POSITION;

    if src.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_i32_le() {
        UNMAPPED_POSITION => Ok(None),
        n => usize::try_from(n + 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .map(Position::new),
    }
}

fn get_flags<B>(src: &mut B) -> io::Result<sam::record::Flags>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(sam::record::Flags::from(src.get_u16_le()))
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
        let mut buf = Vec::new();
        let mut record = Record::default();
        let block_size = read_record(&mut reader, &mut buf, &mut record)?;

        assert_eq!(block_size, 34);
        assert_eq!(record, Record::default());

        Ok(())
    }

    #[test]
    fn test_decode_record_with_invalid_l_read_name() {
        let data = vec![
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x00, // l_read_name = 0
        ];
        let mut src = &data[..];

        let mut record = Record::default();

        assert!(matches!(
            decode_record(&mut src, &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
