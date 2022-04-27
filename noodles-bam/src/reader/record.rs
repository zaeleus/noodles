//! BAM record field readers.

mod cigar;
pub mod data;
mod quality_scores;
mod read_name;
mod sequence;

use std::{
    io::{self, Read},
    mem,
    num::NonZeroUsize,
};

use byteorder::{LittleEndian, ReadBytesExt};
use bytes::{Buf, BytesMut};
use noodles_core::Position;
use noodles_sam as sam;

use crate::Record;

pub(crate) fn read_record<R>(
    reader: &mut R,
    buf: &mut BytesMut,
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

    decode_record(&buf[..], record)?;

    Ok(block_size)
}

pub(crate) fn decode_record<B>(mut src: B, record: &mut Record) -> io::Result<()>
where
    B: Buf,
{
    use self::{
        cigar::get_cigar, data::get_data, quality_scores::get_quality_scores,
        read_name::get_read_name, sequence::get_sequence,
    };

    *record.reference_sequence_id_mut() = get_reference_sequence_id(&mut src)?;
    *record.position_mut() = get_position(&mut src)?;

    let l_read_name = NonZeroUsize::new(usize::from(src.get_u8()))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid l_read_name"))?;

    *record.mapping_quality_mut() = get_mapping_quality(&mut src)?;

    // Discard bin.
    src.advance(mem::size_of::<u16>());

    let n_cigar_op = usize::from(src.get_u16_le());

    *record.flags_mut() = get_flags(&mut src)?;

    let l_seq = usize::try_from(src.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *record.mate_reference_sequence_id_mut() = get_reference_sequence_id(&mut src)?;
    *record.mate_position_mut() = get_position(&mut src)?;

    *record.template_length_mut() = src.get_i32_le();

    get_read_name(&mut src, record.read_name_mut(), l_read_name)?;
    get_cigar(&mut src, record.cigar_mut(), n_cigar_op)?;
    get_sequence(&mut src, record.sequence_mut(), l_seq)?;
    get_quality_scores(&mut src, record.quality_scores_mut(), l_seq)?;

    get_data(&mut src, record.data_mut())?;

    Ok(())
}

fn get_reference_sequence_id<B>(src: &mut B) -> io::Result<Option<usize>>
where
    B: Buf,
{
    use crate::record::reference_sequence_id::UNMAPPED;

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

fn get_mapping_quality<B>(src: &mut B) -> io::Result<Option<sam::record::MappingQuality>>
where
    B: Buf,
{
    use sam::record::mapping_quality::MISSING;

    if src.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_u8() {
        MISSING => Ok(None),
        n => Ok(sam::record::MappingQuality::new(n)),
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
        let mut buf = BytesMut::new();
        let mut record = Record::default();
        let block_size = read_record(&mut reader, &mut buf, &mut record)?;

        assert_eq!(block_size, 34);
        assert_eq!(record, Record::default());

        Ok(())
    }

    #[test]
    fn test_decode_record_with_invalid_l_read_name() {
        let data = [
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x00, // l_read_name = 0
        ];

        let mut reader = &data[..];
        let mut record = Record::default();

        assert!(matches!(
            decode_record(&mut reader, &mut record),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
