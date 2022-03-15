//! BAM record field readers.

pub mod data;

use std::{
    io::{self, Read},
    mem,
    num::NonZeroUsize,
};

use byteorder::{LittleEndian, ReadBytesExt};
use bytes::Buf;
use noodles_core::Position;
use noodles_sam as sam;

use crate::{record::ReferenceSequenceId, Record};

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

    decode_record(&buf[..], record)?;

    Ok(block_size)
}

pub(crate) fn decode_record<B>(mut buf: B, record: &mut Record) -> io::Result<()>
where
    B: Buf,
{
    *record.reference_sequence_id_mut() = read_reference_sequence_id(&mut buf)?;
    *record.position_mut() = read_position(&mut buf)?;

    let l_read_name = NonZeroUsize::new(usize::from(buf.get_u8()))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid l_read_name"))?;

    *record.mapping_quality_mut() = read_mapping_quality(&mut buf)?;

    // Discard bin.
    buf.advance(mem::size_of::<u16>());

    let n_cigar_op = usize::from(buf.get_u16_le());

    *record.flags_mut() = read_flag(&mut buf)?;

    let l_seq = usize::try_from(buf.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *record.mate_reference_sequence_id_mut() = read_reference_sequence_id(&mut buf)?;
    *record.mate_position_mut() = read_position(&mut buf)?;

    *record.template_length_mut() = buf.get_i32_le();

    read_read_name(&mut buf, record.read_name_mut(), l_read_name)?;
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
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match buf.get_i32_le() {
        UNMAPPED => Ok(None),
        n => usize::try_from(n)
            .map(ReferenceSequenceId::from)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

fn read_position<B>(buf: &mut B) -> io::Result<Option<Position>>
where
    B: Buf,
{
    use crate::record::UNMAPPED_POSITION;

    if buf.remaining() < mem::size_of::<i32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match buf.get_i32_le() {
        UNMAPPED_POSITION => Ok(None),
        n => usize::try_from(n + 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .map(Position::new),
    }
}

fn read_mapping_quality<B>(buf: &mut B) -> io::Result<Option<sam::record::MappingQuality>>
where
    B: Buf,
{
    use sam::record::mapping_quality::MISSING;

    if buf.remaining() < mem::size_of::<u8>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
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
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(sam::record::Flags::from(buf.get_u16_le()))
}

fn read_read_name<B>(
    buf: &mut B,
    read_name: &mut Option<sam::record::ReadName>,
    l_read_name: NonZeroUsize,
) -> io::Result<()>
where
    B: Buf,
{
    const MISSING: u8 = b'*';

    let len = usize::from(l_read_name);

    if buf.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    *read_name = if len == 2 && buf.chunk()[0] == MISSING {
        buf.advance(2);
        None
    } else {
        let mut read_name_buf = read_name.take().map(Vec::from).unwrap_or_default();

        // SAFETY: len is guaranteed to be > 0.
        read_name_buf.resize(len - 1, Default::default());
        buf.copy_to_slice(&mut read_name_buf);

        // Discard the NUL terminator.
        buf.advance(1);

        sam::record::ReadName::try_from(read_name_buf)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
    };

    Ok(())
}

fn read_cigar<B>(buf: &mut B, cigar: &mut sam::record::Cigar, n_cigar_op: usize) -> io::Result<()>
where
    B: Buf,
{
    use sam::record::cigar::{op::Kind, Op};

    fn decode_cigar_op(n: u32) -> io::Result<Op> {
        let kind = match n & 0x0f {
            0 => Kind::Match,
            1 => Kind::Insertion,
            2 => Kind::Deletion,
            3 => Kind::Skip,
            4 => Kind::SoftClip,
            5 => Kind::HardClip,
            6 => Kind::Pad,
            7 => Kind::SequenceMatch,
            8 => Kind::SequenceMismatch,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid CIGAR op kind",
                ))
            }
        };

        let len =
            usize::try_from(n >> 4).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        Ok(Op::new(kind, len))
    }

    if buf.remaining() < mem::size_of::<u32>() * n_cigar_op {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    cigar.clear();

    for _ in 0..n_cigar_op {
        let op = decode_cigar_op(buf.get_u32_le())?;
        cigar.as_mut().push(op);
    }

    Ok(())
}

fn read_seq<B>(buf: &mut B, sequence: &mut sam::record::Sequence, l_seq: usize) -> io::Result<()>
where
    B: Buf,
{
    use sam::record::sequence::Base;

    fn decode_bases(n: u8) -> [Base; 2] {
        static BASES: [Base; 16] = [
            Base::Eq,
            Base::A,
            Base::C,
            Base::M,
            Base::G,
            Base::R,
            Base::S,
            Base::V,
            Base::T,
            Base::W,
            Base::Y,
            Base::H,
            Base::K,
            Base::D,
            Base::B,
            Base::N,
        ];

        let l = n >> 4;
        let r = n & 0x0f;

        [BASES[usize::from(l)], BASES[usize::from(r)]]
    }

    let seq_len = (l_seq + 1) / 2;

    if buf.remaining() < seq_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let seq = buf.take(seq_len);
    let bases = seq
        .chunk()
        .iter()
        .copied()
        .flat_map(decode_bases)
        .take(l_seq);

    sequence.clear();

    for base in bases {
        sequence.push(base);
    }

    buf.advance(seq_len);

    Ok(())
}

fn read_qual<B>(
    buf: &mut B,
    quality_scores: &mut sam::record::QualityScores,
    l_seq: usize,
) -> io::Result<()>
where
    B: Buf,
{
    use sam::record::quality_scores::Score;

    if buf.remaining() < l_seq {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    quality_scores.clear();

    let qual = buf.take(l_seq);

    if !is_missing_quality_scores(qual.chunk()) {
        for &b in qual.chunk() {
            let score =
                Score::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            quality_scores.push(score);
        }
    }

    buf.advance(l_seq);

    Ok(())
}

fn is_missing_quality_scores(buf: &[u8]) -> bool {
    use crate::writer::alignment_record::NULL_QUALITY_SCORE;

    buf.iter().all(|&b| b == NULL_QUALITY_SCORE)
}

fn read_data<B>(buf: &mut B, data: &mut sam::record::Data) -> io::Result<()>
where
    B: Buf,
{
    use self::data::get_field;

    data.clear();

    while let Some(field) = get_field(buf)? {
        data.insert(field);
    }

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

    #[test]
    fn test_read_qual_with_sequence_and_no_quality_scores() -> io::Result<()> {
        let data = [0xff, 0xff, 0xff, 0xff];
        let mut buf = &data[..];

        let mut quality_scores = sam::record::QualityScores::default();
        read_qual(&mut buf, &mut quality_scores, 4)?;

        assert!(quality_scores.is_empty());

        Ok(())
    }
}
