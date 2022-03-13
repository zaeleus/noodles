use std::{
    ffi::CString,
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    header::ReferenceSequences,
    record::{sequence::Base, Data},
    AlignmentRecord,
};

// § 4.2 The BAM format (2021-06-03)
//
// ref_id (4) + pos (4) + l_read_name (1) + mapq (1) + bin (2) + n_cigar_op (2) + flag (2) + l_seq
// (4) + next_ref_id (4) + next_pos (4) + tlen (4)
const BLOCK_HEADER_SIZE: u32 = 32;

// § 4.2.3 SEQ and QUAL encoding (2021-06-03)
pub(crate) const NULL_QUALITY_SCORE: u8 = 255;

pub fn write_sam_record<W>(
    writer: &mut W,
    reference_sequences: &ReferenceSequences,
    record: &sam::Record,
) -> io::Result<()>
where
    W: Write,
{
    let name = record.read_name().map(|name| name.as_ref()).unwrap_or("*");
    let c_read_name =
        CString::new(name).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let read_name = c_read_name.as_bytes_with_nul();
    let l_read_name = u8::try_from(read_name.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let n_cigar_op = u16::try_from(record.cigar().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let data_len = calculate_data_len(record.data()).and_then(|n| {
        u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;

    let block_size = BLOCK_HEADER_SIZE
        + u32::from(l_read_name)
        + (4 * u32::from(n_cigar_op))
        + ((l_seq + 1) / 2)
        + l_seq
        + data_len;
    writer.write_u32::<LittleEndian>(block_size)?;

    // ref_id
    write_reference_sequence_id(
        writer,
        reference_sequences,
        record.reference_sequence_name(),
    )?;

    // pos
    write_position(writer, record.alignment_start())?;

    writer.write_u8(l_read_name)?;

    // mapq
    write_mapping_quality(writer, record.mapping_quality())?;

    // bin
    write_bin(writer, record.alignment_start(), record.alignment_end())?;

    writer.write_u16::<LittleEndian>(n_cigar_op)?;

    // flag
    write_flags(writer, record.flags())?;

    writer.write_u32::<LittleEndian>(l_seq)?;

    // next_ref_id
    write_reference_sequence_id(
        writer,
        reference_sequences,
        record.mate_reference_sequence_name(),
    )?;

    // next_pos
    write_position(writer, record.mate_alignment_start())?;

    // tlen
    write_template_length(writer, record.template_length())?;

    writer.write_all(read_name)?;

    write_cigar(writer, record.cigar())?;

    // § 4.2.3 SEQ and QUAL encoding (2021-06-03)
    let sequence = record.sequence();
    let quality_scores = record.quality_scores();

    write_sequence(writer, sequence)?;

    if sequence.len() == quality_scores.len() {
        write_quality_scores(writer, quality_scores)?;
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

    write_data(writer, record.data())?;

    Ok(())
}

fn write_reference_sequence_id<W>(
    writer: &mut W,
    reference_sequences: &ReferenceSequences,
    reference_sequence_name: Option<&sam::record::ReferenceSequenceName>,
) -> io::Result<()>
where
    W: Write,
{
    use crate::record::reference_sequence_id;

    let id = match reference_sequence_name {
        Some(name) => reference_sequences
            .get_index_of(name.as_str())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence name: {}", name),
                )
            })
            .and_then(|i| {
                i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?,
        None => reference_sequence_id::UNMAPPED,
    };

    writer.write_i32::<LittleEndian>(id)
}

fn write_position<W>(writer: &mut W, position: Option<Position>) -> io::Result<()>
where
    W: Write,
{
    use crate::record::UNMAPPED_POSITION;

    let pos = if let Some(position) = position {
        i32::try_from(usize::from(position) - 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
    } else {
        UNMAPPED_POSITION
    };

    writer.write_i32::<LittleEndian>(pos)
}

fn write_mapping_quality<W>(
    writer: &mut W,
    mapping_quality: Option<sam::record::MappingQuality>,
) -> io::Result<()>
where
    W: Write,
{
    use sam::record::mapping_quality::MISSING;
    let mapq = mapping_quality.map(u8::from).unwrap_or(MISSING);
    writer.write_u8(mapq)
}

fn write_bin<W>(
    writer: &mut W,
    alignment_start: Option<Position>,
    alignment_end: Option<Position>,
) -> io::Result<()>
where
    W: Write,
{
    use super::record::{region_to_bin, UNMAPPED_BIN};

    let bin = match (alignment_start, alignment_end) {
        (Some(start), Some(end)) => region_to_bin(start, end)?,
        _ => UNMAPPED_BIN,
    };

    writer.write_u16::<LittleEndian>(bin)
}

fn write_flags<W>(writer: &mut W, flags: sam::record::Flags) -> io::Result<()>
where
    W: Write,
{
    let flag = u16::from(flags);
    writer.write_u16::<LittleEndian>(flag)
}

fn write_template_length<W>(writer: &mut W, template_length: i32) -> io::Result<()>
where
    W: Write,
{
    writer.write_i32::<LittleEndian>(template_length)
}

fn write_cigar<W>(writer: &mut W, cigar: &sam::record::Cigar) -> io::Result<()>
where
    W: Write,
{
    use super::record::encode_cigar_op;

    for &op in cigar.as_ref() {
        let n = encode_cigar_op(op)?;
        writer.write_u32::<LittleEndian>(n)?;
    }

    Ok(())
}

fn write_sequence<W>(writer: &mut W, sequence: &sam::record::Sequence) -> io::Result<()>
where
    W: Write,
{
    use super::record::encode_base;

    for chunk in sequence.as_ref().chunks(2) {
        let l = chunk[0];
        let r = chunk.get(1).copied().unwrap_or(Base::Eq);
        let b = encode_base(l) << 4 | encode_base(r);
        writer.write_u8(b)?;
    }

    Ok(())
}

fn write_quality_scores<W>(
    writer: &mut W,
    quality_scores: &sam::record::QualityScores,
) -> io::Result<()>
where
    W: Write,
{
    for &score in quality_scores.iter() {
        writer.write_u8(u8::from(score))?;
    }

    Ok(())
}

pub(crate) fn calculate_data_len(data: &Data) -> io::Result<usize> {
    use noodles_sam::record::data::field::Value;

    let mut len = 0;

    for field in data.values() {
        // tag
        len += 2;
        // val_type
        len += 1;

        let value = field.value();

        if value.subtype().is_some() {
            // subtype
            len += 1;
            // count
            len += mem::size_of::<u32>();
        }

        len += match value {
            Value::Char(_) => mem::size_of::<u8>(),
            Value::Int8(_) => mem::size_of::<i8>(),
            Value::UInt8(_) => mem::size_of::<u8>(),
            Value::Int16(_) => mem::size_of::<i16>(),
            Value::UInt16(_) => mem::size_of::<u16>(),
            Value::Int32(_) => mem::size_of::<i32>(),
            Value::UInt32(_) => mem::size_of::<u32>(),
            Value::Float(_) => mem::size_of::<f32>(),
            Value::String(s) | Value::Hex(s) => s.as_bytes().len() + 1,
            Value::Int8Array(values) => values.len(),
            Value::UInt8Array(values) => values.len(),
            Value::Int16Array(values) => mem::size_of::<i16>() * values.len(),
            Value::UInt16Array(values) => mem::size_of::<u16>() * values.len(),
            Value::Int32Array(values) => mem::size_of::<i32>() * values.len(),
            Value::UInt32Array(values) => mem::size_of::<u32>() * values.len(),
            Value::FloatArray(values) => mem::size_of::<f32>() * values.len(),
        }
    }

    Ok(len)
}

fn write_data<W>(writer: &mut W, data: &Data) -> io::Result<()>
where
    W: Write,
{
    use crate::writer::record::data::write_field;

    for field in data.values() {
        write_field(writer, field)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_sequence_id() -> Result<(), Box<dyn std::error::Error>> {
        use sam::header::{reference_sequence, ReferenceSequence};

        let reference_sequences = [("sq0".parse()?, 8), ("sq1".parse()?, 13)]
            .into_iter()
            .map(|(name, len): (reference_sequence::Name, i32)| {
                let sn = name.to_string();
                ReferenceSequence::new(name, len).map(|rs| (sn, rs))
            })
            .collect::<Result<_, _>>()?;

        let mut buf = Vec::new();

        buf.clear();
        let reference_sequence_name = "sq0".parse()?;
        write_reference_sequence_id(
            &mut buf,
            &reference_sequences,
            Some(&reference_sequence_name),
        )?;
        assert_eq!(buf, [0x00, 0x00, 0x00, 0x00]);

        buf.clear();
        write_reference_sequence_id(&mut buf, &reference_sequences, None)?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff]);

        buf.clear();
        let reference_sequence_name = "sq2".parse()?;
        assert!(matches!(
            write_reference_sequence_id(&mut buf, &reference_sequences, Some(&reference_sequence_name)),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }
}
