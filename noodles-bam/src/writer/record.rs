use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use super::sam_record::NULL_QUALITY_SCORE;
use crate::{record::ReferenceSequenceId, Record};

pub(super) fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    let block_size = u32::try_from(record.block_size())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(block_size)?;

    // ref_id
    write_reference_sequence_id(writer, record.reference_sequence_id())?;

    // pos
    write_position(writer, record.pos)?;

    let l_read_name = u8::try_from(record.read_name.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u8(l_read_name)?;

    let mapq = u8::from(record.mapping_quality());
    writer.write_u8(mapq)?;

    writer.write_u16::<LittleEndian>(record.bin())?;

    let n_cigar_op = u16::try_from(record.cigar().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(n_cigar_op)?;

    let flag = u16::from(record.flags());
    writer.write_u16::<LittleEndian>(flag)?;

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(l_seq)?;

    // next_ref_id
    write_reference_sequence_id(writer, record.mate_reference_sequence_id())?;

    // next_pos
    write_position(writer, record.next_pos)?;

    writer.write_i32::<LittleEndian>(record.template_length())?;

    writer.write_all(&record.read_name)?;

    for &raw_op in record.cigar().as_ref().iter() {
        writer.write_u32::<LittleEndian>(raw_op)?;
    }

    let sequence = record.sequence();
    let quality_scores = record.quality_scores();

    writer.write_all(sequence.as_ref())?;

    if sequence.len() == quality_scores.len() {
        writer.write_all(record.quality_scores().as_ref())?;
    } else if quality_scores.is_empty() {
        for _ in 0..sequence.len() {
            writer.write_u8(NULL_QUALITY_SCORE)?;
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

    writer.write_all(record.data().as_ref())?;

    Ok(())
}

fn write_reference_sequence_id<W>(
    writer: &mut W,
    reference_sequence_id: Option<ReferenceSequenceId>,
) -> io::Result<()>
where
    W: Write,
{
    use crate::record::reference_sequence_id::UNMAPPED;

    let ref_id = reference_sequence_id.map(i32::from).unwrap_or(UNMAPPED);
    writer.write_i32::<LittleEndian>(ref_id)
}

// pos is 0-based.
fn write_position<W>(writer: &mut W, pos: i32) -> io::Result<()>
where
    W: Write,
{
    writer.write_i32::<LittleEndian>(pos)
}

#[cfg(test)]
mod tests {
    use noodles_sam as sam;

    use super::*;

    #[test]
    fn test_write_record_with_default_fields() -> Result<(), Box<dyn std::error::Error>> {
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

    #[test]
    fn test_write_record_with_all_fields() -> Result<(), Box<dyn std::error::Error>> {
        use sam::record::{
            cigar::op::Kind, data::field::Tag, quality_scores::Score, Flags, MappingQuality,
        };

        use crate::record::{
            cigar::Op,
            data::{field::Value, Field},
            sequence::Base,
            ReferenceSequenceId,
        };

        let mut record = Record::default();

        *record.reference_sequence_id_mut() = ReferenceSequenceId::try_from(1).map(Some)?;
        record.pos = 8; // 0-based
        *record.mapping_quality_mut() = MappingQuality::from(13);
        *record.bin_mut() = 6765;
        *record.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        *record.mate_reference_sequence_id_mut() = record.reference_sequence_id();
        record.next_pos = 21; // 0-based
        *record.template_length_mut() = 144;

        record.read_name.clear();
        record.read_name.extend_from_slice(b"r0\x00");

        let cigar = record.cigar_mut();
        cigar.push(Op::new(Kind::Match, 36)?);
        cigar.push(Op::new(Kind::SoftClip, 8)?);

        let sequence = record.sequence_mut();
        sequence.push(Base::A);
        sequence.push(Base::C);
        sequence.push(Base::G);
        sequence.push(Base::T);

        let quality_scores = record.quality_scores_mut();
        quality_scores.push(Score::try_from('N')?);
        quality_scores.push(Score::try_from('D')?);
        quality_scores.push(Score::try_from('L')?);
        quality_scores.push(Score::try_from('S')?);

        let data = record.data_mut();
        data.insert(Field::new(Tag::AlignmentHitCount, Value::UInt8(1)));

        let mut buf = Vec::new();
        write_record(&mut buf, &record)?;

        let expected = [
            0x35, 0x00, 0x00, 0x00, // block_size = 53
            0x01, 0x00, 0x00, 0x00, // ref_id = 1
            0x08, 0x00, 0x00, 0x00, // pos = 8
            0x03, // l_read_name = 3
            0x0d, // mapq = 13
            0x6d, 0x1a, // bin = 6765
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
}
