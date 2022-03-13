pub(crate) mod data;

pub(crate) use self::data::put_data;

use std::io;

use bytes::BufMut;
use noodles_core::Position;
use noodles_sam::{self as sam, record::sequence::Base, AlignmentRecord};

use super::sam_record::NULL_QUALITY_SCORE;
use crate::{record::ReferenceSequenceId, Record};

// § 4.2.1 "BIN field calculation" (2021-06-03): "Note unmapped reads with `POS` 0 (which
// becomes -1 in BAM) therefore use `reg2bin(-1, 0)` which is computed as 4680."
pub(crate) const UNMAPPED_BIN: u16 = 4680;

pub(crate) fn encode_record<B>(dst: &mut B, record: &Record) -> io::Result<()>
where
    B: BufMut,
{
    // ref_id
    put_reference_sequence_id(dst, record.reference_sequence_id())?;

    // pos
    put_position(dst, record.alignment_start())?;

    put_l_read_name(dst, record.read_name())?;

    // mapq
    put_mapping_quality(dst, record.mapping_quality());

    // bin
    put_bin(dst, record.alignment_start(), record.alignment_end())?;

    let n_cigar_op = u16::try_from(record.cigar().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u16_le(n_cigar_op);

    // flag
    put_flags(dst, record.flags());

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u32_le(l_seq);

    // next_ref_id
    put_reference_sequence_id(dst, record.mate_reference_sequence_id())?;

    // next_pos
    put_position(dst, record.mate_alignment_start())?;

    // tlen
    put_template_length(dst, record.template_length());

    put_read_name(dst, record.read_name());

    put_cigar(dst, record.cigar())?;

    let sequence = record.sequence();
    let quality_scores = record.quality_scores();

    // seq
    put_sequence(dst, sequence);

    if sequence.len() == quality_scores.len() {
        put_quality_scores(dst, quality_scores);
    } else if quality_scores.is_empty() {
        for _ in 0..sequence.len() {
            dst.put_u8(NULL_QUALITY_SCORE);
        }
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "quality scores length mismatch: expected {}, got {}",
                sequence.len(),
                quality_scores.len()
            ),
        ));
    }

    put_data(dst, record.data())?;

    Ok(())
}

fn put_reference_sequence_id<B>(
    dst: &mut B,
    reference_sequence_id: Option<ReferenceSequenceId>,
) -> io::Result<()>
where
    B: BufMut,
{
    use crate::record::reference_sequence_id::UNMAPPED;

    let ref_id = if let Some(id) = reference_sequence_id {
        i32::try_from(usize::from(id))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
    } else {
        UNMAPPED
    };

    dst.put_i32_le(ref_id);

    Ok(())
}

pub(super) fn put_position<B>(dst: &mut B, position: Option<Position>) -> io::Result<()>
where
    B: BufMut,
{
    use crate::record::UNMAPPED_POSITION;

    let pos = if let Some(position) = position {
        i32::try_from(usize::from(position) - 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
    } else {
        UNMAPPED_POSITION
    };

    dst.put_i32_le(pos);

    Ok(())
}

pub(super) fn put_l_read_name<B>(
    dst: &mut B,
    read_name: Option<&sam::record::ReadName>,
) -> io::Result<()>
where
    B: BufMut,
{
    use std::mem;

    let mut read_name_len = read_name
        .map(|name| name.len())
        .unwrap_or(sam::record::read_name::MISSING.len());

    // + NUL terminator
    read_name_len += mem::size_of::<u8>();

    let l_read_name =
        u8::try_from(read_name_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    dst.put_u8(l_read_name);

    Ok(())
}

pub(super) fn put_mapping_quality<B>(
    dst: &mut B,
    mapping_quality: Option<sam::record::MappingQuality>,
) where
    B: BufMut,
{
    use sam::record::mapping_quality::MISSING;
    let mapq = mapping_quality.map(u8::from).unwrap_or(MISSING);
    dst.put_u8(mapq);
}

pub(super) fn put_bin<B>(
    dst: &mut B,
    alignment_start: Option<Position>,
    alignment_end: Option<Position>,
) -> io::Result<()>
where
    B: BufMut,
{
    let bin = match (alignment_start, alignment_end) {
        (Some(start), Some(end)) => region_to_bin(start, end)?,
        _ => UNMAPPED_BIN,
    };

    dst.put_u16_le(bin);

    Ok(())
}

pub(super) fn put_flags<B>(dst: &mut B, flags: sam::record::Flags)
where
    B: BufMut,
{
    let flag = u16::from(flags);
    dst.put_u16_le(flag);
}

pub(super) fn put_template_length<B>(dst: &mut B, template_length: i32)
where
    B: BufMut,
{
    dst.put_i32_le(template_length);
}

pub(super) fn put_read_name<B>(dst: &mut B, read_name: Option<&sam::record::ReadName>)
where
    B: BufMut,
{
    use sam::record::read_name::MISSING;

    const NUL: u8 = 0x00;

    if let Some(read_name) = read_name {
        dst.put(read_name.as_ref());
    } else {
        dst.put(MISSING);
    }

    dst.put_u8(NUL);
}

pub(super) fn put_cigar<B>(dst: &mut B, cigar: &sam::record::Cigar) -> io::Result<()>
where
    B: BufMut,
{
    for &op in cigar.as_ref() {
        let n = encode_cigar_op(op)?;
        dst.put_u32_le(n);
    }

    Ok(())
}

pub(crate) fn encode_cigar_op(op: sam::record::cigar::Op) -> io::Result<u32> {
    const MAX_LENGTH: u32 = (1 << 28) - 1;

    let len =
        u32::try_from(op.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if len <= MAX_LENGTH {
        let k = op.kind() as u32;
        Ok(len << 4 | k)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CIGAR op length",
        ))
    }
}

pub(super) fn put_sequence<B>(dst: &mut B, sequence: &sam::record::Sequence)
where
    B: BufMut,
{
    for chunk in sequence.as_ref().chunks(2) {
        let l = chunk[0];
        let r = chunk.get(1).copied().unwrap_or(Base::Eq);
        let b = encode_base(l) << 4 | encode_base(r);
        dst.put_u8(b);
    }
}

pub(crate) fn encode_base(base: Base) -> u8 {
    match base {
        Base::Eq => 0,
        Base::A => 1,
        Base::C => 2,
        Base::M => 3,
        Base::G => 4,
        Base::R => 5,
        Base::S => 6,
        Base::V => 7,
        Base::T => 8,
        Base::W => 9,
        Base::Y => 10,
        Base::H => 11,
        Base::K => 12,
        Base::D => 13,
        Base::B => 14,
        // § 4.2.3 SEQ and QUAL encoding (2021-06-03): "The case-insensitive base codes ... are
        // mapped to [0, 15] respectively with all other characters mapping to 'N' (value 15)".
        _ => 15,
    }
}

pub(super) fn put_quality_scores<B>(dst: &mut B, quality_scores: &sam::record::QualityScores)
where
    B: BufMut,
{
    for &score in quality_scores.iter() {
        dst.put_u8(u8::from(score));
    }
}

// § 5.3 "C source code for computing bin number and overlapping bins" (2021-06-03)
#[allow(clippy::eq_op)]
pub(crate) fn region_to_bin(alignment_start: Position, alignment_end: Position) -> io::Result<u16> {
    let start = usize::from(alignment_start) - 1;
    let end = usize::from(alignment_end) - 1;

    let bin = if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (start >> 26)
    } else {
        0
    };

    u16::try_from(bin).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record_with_default_fields() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let record = Record::default();
        encode_record(&mut buf, &record)?;

        let expected = [
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

    #[test]
    fn test_write_record_with_all_fields() -> Result<(), Box<dyn std::error::Error>> {
        use sam::record::{
            data::{
                field::{Tag, Value},
                Field,
            },
            Data, Flags, MappingQuality,
        };

        use crate::record::ReferenceSequenceId;

        let reference_sequence_id = ReferenceSequenceId::from(1);

        let record = Record::builder()
            .set_reference_sequence_id(reference_sequence_id)
            .set_position(Position::try_from(9)?)
            .set_mapping_quality(MappingQuality::try_from(13)?)
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_mate_reference_sequence_id(reference_sequence_id)
            .set_mate_position(Position::try_from(22)?)
            .set_template_length(144)
            .set_read_name("r0".parse()?)
            .set_cigar("36M8S".parse()?)
            .set_sequence("ACGT".parse()?)
            .set_quality_scores("NDLS".parse()?)
            .set_data(Data::try_from(vec![Field::new(
                Tag::AlignmentHitCount,
                Value::UInt8(1),
            )])?)
            .build()?;

        let mut buf = Vec::new();
        encode_record(&mut buf, &record)?;

        let expected = [
            0x01, 0x00, 0x00, 0x00, // ref_id = 1
            0x08, 0x00, 0x00, 0x00, // pos = 8
            0x03, // l_read_name = 3
            0x0d, // mapq = 13
            0x49, 0x12, // bin = 4681
            0x02, 0x00, // n_cigar_op = 2
            0x41, 0x00, // flag = 65
            0x04, 0x00, 0x00, 0x00, // l_seq = 4
            0x01, 0x00, 0x00, 0x00, // next_ref_id = 1
            0x15, 0x00, 0x00, 0x00, // next_pos = 21
            0x90, 0x00, 0x00, 0x00, // tlen = 144
            b'r', b'0', 0x00, // read_name = "r0\x00"
            0x40, 0x02, 0x00, 0x00, // cigar[0] = 36M
            0x84, 0x00, 0x00, 0x00, // cigar[1] = 8S
            0x12, 0x48, // seq = ACGT
            0x2d, 0x23, 0x2b, 0x32, // qual = NDLS
            b'N', b'H', b'C', 0x01, // data[0] = NH:i:1
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_region_to_bin() -> Result<(), Box<dyn std::error::Error>> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        assert_eq!(region_to_bin(start, end)?, 4681);

        let start = Position::try_from(63245986)?;
        let end = Position::try_from(63245986)?;
        assert_eq!(region_to_bin(start, end)?, 8541);

        Ok(())
    }
}
