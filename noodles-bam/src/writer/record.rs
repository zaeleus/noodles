use std::{
    convert::TryFrom,
    ffi::CString,
    io::{self, Write},
    mem, num,
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_sam::{
    self as sam,
    header::ReferenceSequences,
    record::{Cigar, Data, QualityScores, Sequence},
};

use crate::record::sequence::Base;

// § 4.2 The BAM format (2021-06-03)
//
// ref_id (4) + pos (4) + l_read_name (1) + mapq (1) + bin (2) + n_cigar_op (2) + flag (2) + l_seq
// (4) + next_ref_id (4) + next_pos (4) + tlen (4)
const BLOCK_HEADER_SIZE: u32 = 32;

// § 4.2.3 SEQ and QUAL encoding (2021-06-03)
const NULL_QUALITY_SCORE: u8 = 255;

pub fn write_sam_record<W>(
    writer: &mut W,
    reference_sequences: &ReferenceSequences,
    record: &sam::Record,
) -> io::Result<()>
where
    W: Write,
{
    let name = record.read_name().map(|name| name.as_str()).unwrap_or("*");
    let c_read_name =
        CString::new(name).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let read_name = c_read_name.as_bytes_with_nul();
    let l_read_name = u8::try_from(read_name.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let n_cigar_op = u16::try_from(record.cigar().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let data_len = calculate_data_len(record.data()).and_then(|len| {
        u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
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
    write_position(writer, record.position())?;

    writer.write_u8(l_read_name)?;

    // mapq
    write_mapping_quality(writer, record.mapping_quality())?;

    // bin
    write_bin(writer, record)?;

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
    write_position(writer, record.mate_position())?;

    // tlen
    write_template_length(writer, record.template_length())?;

    writer.write_all(read_name)?;

    write_cigar(writer, record.cigar())?;

    // § 4.2.3 SEQ and QUAL encoding (2021-06-03)
    let sequence = record.sequence();
    let quality_scores = record.quality_scores();

    if !sequence.is_empty() {
        write_seq(writer, sequence)?;

        if sequence.len() == quality_scores.len() {
            write_qual(writer, quality_scores)?;
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

fn write_position<W>(writer: &mut W, position: Option<sam::record::Position>) -> io::Result<()>
where
    W: Write,
{
    use crate::record::UNMAPPED_POSITION;

    let pos = position
        .map(|p| i32::from(p) - 1)
        .unwrap_or(UNMAPPED_POSITION);

    writer.write_i32::<LittleEndian>(pos)
}

fn write_mapping_quality<W>(
    writer: &mut W,
    mapping_quality: sam::record::MappingQuality,
) -> io::Result<()>
where
    W: Write,
{
    let mapq = u8::from(mapping_quality);
    writer.write_u8(mapq)
}

fn write_bin<W>(writer: &mut W, record: &sam::Record) -> io::Result<()>
where
    W: Write,
{
    // § 4.2.1 "BIN field calculation" (2021-06-03): "Note unmapped reads with `POS` 0 (which
    // becomes -1 in BAM) therefore use `reg2bin(-1, 0)` which is computed as 4680."
    const UNMAPPED_BIN: u16 = 4680;

    let bin = record
        .position()
        .map(|p| i32::from(p) - 1)
        .map(|start| {
            let reference_len = record.cigar().reference_len() as i32;
            let end = start + reference_len;
            region_to_bin(start, end) as u16
        })
        .unwrap_or(UNMAPPED_BIN);

    writer.write_u16::<LittleEndian>(bin)
}

fn write_flags<W>(writer: &mut W, flags: sam::record::Flags) -> io::Result<()>
where
    W: Write,
{
    let flag = u16::from(flags);
    writer.write_u16::<LittleEndian>(flag)
}

fn write_template_length<W>(writer: &mut W, tlen: i32) -> io::Result<()>
where
    W: Write,
{
    writer.write_i32::<LittleEndian>(tlen)
}

fn write_cigar<W>(writer: &mut W, cigar: &Cigar) -> io::Result<()>
where
    W: Write,
{
    for op in cigar.iter() {
        let len = op.len();
        let kind = op.kind() as u32;
        let value = len << 4 | kind;
        writer.write_u32::<LittleEndian>(value)?;
    }

    Ok(())
}

fn write_seq<W>(writer: &mut W, sequence: &Sequence) -> io::Result<()>
where
    W: Write,
{
    for chunk in sequence.chunks(2) {
        let l = Base::from(chunk[0]);
        let r = chunk.get(1).copied().map(Base::from).unwrap_or(Base::Eq);
        let value = u8::from(l) << 4 | u8::from(r);
        writer.write_u8(value)?;
    }

    Ok(())
}

fn write_qual<W>(writer: &mut W, quality_scores: &QualityScores) -> io::Result<()>
where
    W: Write,
{
    for score in quality_scores.iter() {
        let value = u8::from(*score);
        writer.write_u8(value)?;
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

        match value {
            Value::Char(_) => {
                len += mem::size_of::<u8>();
            }
            Value::Int(n) => {
                if *n >= 0 {
                    if *n <= i64::from(u8::MAX) {
                        len += mem::size_of::<u8>();
                    } else if *n <= i64::from(u16::MAX) {
                        len += mem::size_of::<u16>();
                    } else if *n <= i64::from(u32::MAX) {
                        len += mem::size_of::<u32>();
                    } else {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("invalid integer value: {}", n),
                        ));
                    }
                } else if *n >= i64::from(i8::MIN) {
                    len += mem::size_of::<i8>();
                } else if *n >= i64::from(i16::MIN) {
                    len += mem::size_of::<i16>();
                } else if *n >= i64::from(i32::MIN) {
                    len += mem::size_of::<i32>();
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid integer value: {}", n),
                    ));
                }
            }
            Value::Float(_) => {
                len += mem::size_of::<f32>();
            }
            Value::String(s) | Value::Hex(s) => {
                len += s.as_bytes().len() + 1;
            }
            Value::Int8Array(values) => {
                len += values.len();
            }
            Value::UInt8Array(values) => {
                len += values.len();
            }
            Value::Int16Array(values) => {
                len += mem::size_of::<i16>() * values.len();
            }
            Value::UInt16Array(values) => {
                len += mem::size_of::<u16>() * values.len();
            }
            Value::Int32Array(values) => {
                len += mem::size_of::<i32>() * values.len();
            }
            Value::UInt32Array(values) => {
                len += mem::size_of::<u32>() * values.len();
            }
            Value::FloatArray(values) => {
                len += mem::size_of::<f32>() * values.len();
            }
        }
    }

    Ok(len)
}

fn write_data<W>(writer: &mut W, data: &Data) -> io::Result<()>
where
    W: Write,
{
    use noodles_sam::record::data::field::Value;

    for field in data.values() {
        let tag = field.tag();
        writer.write_all(tag.as_ref())?;

        let value = field.value();

        if let Value::Int(n) = value {
            write_data_field_int_value(writer, *n)?;
        } else {
            write_data_field_value_type(writer, value)?;
            write_data_field_value(writer, value)?;
        }
    }

    Ok(())
}

fn write_data_field_int_value<W>(writer: &mut W, n: i64) -> io::Result<()>
where
    W: Write,
{
    use crate::record::data::field::value::Type;

    if n >= 0 {
        if n <= i64::from(u8::MAX) {
            writer.write_u8(u8::from(Type::UInt8))?;
            return writer.write_u8(n as u8);
        } else if n <= i64::from(u16::MAX) {
            writer.write_u8(u8::from(Type::UInt16))?;
            return writer.write_u16::<LittleEndian>(n as u16);
        } else if n <= i64::from(u32::MAX) {
            writer.write_u8(u8::from(Type::UInt32))?;
            return writer.write_u32::<LittleEndian>(n as u32);
        }
    } else if n >= i64::from(i8::MIN) {
        writer.write_u8(u8::from(Type::Int8))?;
        return writer.write_i8(n as i8);
    } else if n >= i64::from(i16::MIN) {
        writer.write_u8(u8::from(Type::Int16))?;
        return writer.write_i16::<LittleEndian>(n as i16);
    } else if n >= i64::from(i32::MIN) {
        writer.write_u8(u8::from(Type::Int32))?;
        return writer.write_i32::<LittleEndian>(n as i32);
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("invalid integer value: {}", n),
    ))
}

fn write_data_field_value_type<W>(
    writer: &mut W,
    value: &sam::record::data::field::Value,
) -> io::Result<()>
where
    W: Write,
{
    let val_type = char::from(value.ty()) as u8;
    writer.write_u8(val_type)?;

    if let Some(subtype) = value.subtype() {
        let val_subtype = char::from(subtype) as u8;
        writer.write_u8(val_subtype)?;
    }

    Ok(())
}

fn write_data_field_value<W>(
    writer: &mut W,
    value: &sam::record::data::field::Value,
) -> io::Result<()>
where
    W: Write,
{
    use sam::record::data::field::Value;

    fn invalid_array_len(e: num::TryFromIntError) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidInput, e)
    }

    match value {
        Value::Char(c) => writer.write_u8(*c as u8)?,
        Value::Int(_) => {
            // Integers are handled by `write_data_field_int_value`.
            unreachable!();
        }
        Value::Float(n) => writer.write_f32::<LittleEndian>(*n)?,
        Value::String(s) | Value::Hex(s) => {
            let c_str = CString::new(s.as_bytes())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            writer.write_all(c_str.as_bytes_with_nul())?;
        }
        Value::Int8Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32::<LittleEndian>(len)?;

            for &n in values {
                writer.write_i8(n)?;
            }
        }
        Value::UInt8Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32::<LittleEndian>(len)?;

            for &n in values {
                writer.write_u8(n)?;
            }
        }
        Value::Int16Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32::<LittleEndian>(len)?;

            for &n in values {
                writer.write_i16::<LittleEndian>(n)?;
            }
        }
        Value::UInt16Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32::<LittleEndian>(len)?;

            for &n in values {
                writer.write_u16::<LittleEndian>(n)?;
            }
        }
        Value::Int32Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32::<LittleEndian>(len)?;

            for &n in values {
                writer.write_i32::<LittleEndian>(n)?;
            }
        }
        Value::UInt32Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32::<LittleEndian>(len)?;

            for &n in values {
                writer.write_u32::<LittleEndian>(n)?;
            }
        }
        Value::FloatArray(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32::<LittleEndian>(len)?;

            for &n in values {
                writer.write_f32::<LittleEndian>(n)?;
            }
        }
    }

    Ok(())
}

// § 5.3 C source code for computing bin number and overlapping bins (2021-06-03)
// 0-based, [start, end)
#[allow(clippy::eq_op)]
pub(crate) fn region_to_bin(start: i32, mut end: i32) -> i32 {
    end -= 1;

    if start >> 14 == end >> 14 {
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
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_data_field_int_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: i64, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data_field_int_value(buf, n)?;
            assert_eq!(&buf[..], expected, "n = {}", n);
            Ok(())
        }

        let mut buf = Vec::new();

        // i32::MIN - 1
        buf.clear();
        assert!(matches!(
            write_data_field_int_value(&mut buf, -2147483649),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));
        // i32::MIN
        t(&mut buf, -2147483648, &[b'i', 0x00, 0x00, 0x00, 0x80])?;
        // i32::MIN + 1
        t(&mut buf, -2147483647, &[b'i', 0x01, 0x00, 0x00, 0x80])?;

        // i16::MIN - 1
        t(&mut buf, -32769, &[b'i', 0xff, 0x7f, 0xff, 0xff])?;
        // i16::MIN
        t(&mut buf, -32768, &[b's', 0x00, 0x80])?;
        // i16::MIN + 1
        t(&mut buf, -32767, &[b's', 0x01, 0x80])?;

        // i8::MIN - 1
        t(&mut buf, -129, &[b's', 0x7f, 0xff])?;
        // i8::MIN
        t(&mut buf, -128, &[b'c', 0x80])?;
        // i8::MIN + 1
        t(&mut buf, -127, &[b'c', 0x81])?;

        // -1
        t(&mut buf, -1, &[b'c', 0xff])?;
        // 0
        t(&mut buf, 0, &[b'C', 0x00])?;
        // 1
        t(&mut buf, 1, &[b'C', 0x01])?;

        // i8::MAX - 1
        t(&mut buf, 126, &[b'C', 0x7e])?;
        // i8::MAX
        t(&mut buf, 127, &[b'C', 0x7f])?;
        // i8::MAX + 1
        t(&mut buf, 128, &[b'C', 0x80])?;

        // u8::MAX - 1
        t(&mut buf, 254, &[b'C', 0xfe])?;
        // u8::MAX
        t(&mut buf, 255, &[b'C', 0xff])?;
        // u8::MAX + 1
        t(&mut buf, 256, &[b'S', 0x00, 0x01])?;

        // i16::MAX - 1
        t(&mut buf, 32766, &[b'S', 0xfe, 0x7f])?;
        // i16::MAX
        t(&mut buf, 32767, &[b'S', 0xff, 0x7f])?;
        // i16::MAX + 1
        t(&mut buf, 32768, &[b'S', 0x00, 0x80])?;

        // u16::MAX - 1
        t(&mut buf, 65534, &[b'S', 0xfe, 0xff])?;
        // u16::MAX
        t(&mut buf, 65535, &[b'S', 0xff, 0xff])?;
        // u16::MAX + 1
        t(&mut buf, 65536, &[b'I', 0x00, 0x00, 0x01, 0x00])?;

        // i32::MAX - 1
        t(&mut buf, 2147483646, &[b'I', 0xfe, 0xff, 0xff, 0x7f])?;
        // i32::MAX
        t(&mut buf, 2147483647, &[b'I', 0xff, 0xff, 0xff, 0x7f])?;
        // i32::MAX + 1
        t(&mut buf, 2147483648, &[b'I', 0x00, 0x00, 0x00, 0x80])?;

        // u32::MAX - 1
        t(&mut buf, 4294967295, &[b'I', 0xff, 0xff, 0xff, 0xff])?;
        // u32::MAX
        buf.clear();
        assert!(matches!(
            write_data_field_int_value(&mut buf, 4294967296),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_data_field_value_type() -> io::Result<()> {
        use sam::record::data::field::Value;

        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data_field_value_type(buf, value)?;
            assert_eq!(&buf[..], expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Char('n'), &[b'A'])?;
        t(&mut buf, &Value::Int(13), &[b'i'])?;
        t(&mut buf, &Value::Float(8.0), &[b'f'])?;
        t(&mut buf, &Value::String(String::from("ndls")), &[b'Z'])?;
        t(&mut buf, &Value::Hex(String::from("cafe")), &[b'H'])?;
        t(&mut buf, &Value::Int8Array(vec![1, -2]), &[b'B', b'c'])?;
        t(&mut buf, &Value::UInt8Array(vec![3, 5]), &[b'B', b'C'])?;
        t(&mut buf, &Value::Int16Array(vec![8, -13]), &[b'B', b's'])?;
        t(&mut buf, &Value::UInt16Array(vec![21, 34]), &[b'B', b'S'])?;
        t(&mut buf, &Value::Int32Array(vec![55, -89]), &[b'B', b'i'])?;
        t(&mut buf, &Value::UInt32Array(vec![144, 223]), &[b'B', b'I'])?;
        t(&mut buf, &Value::FloatArray(vec![8.0, 13.0]), &[b'B', b'f'])?;

        Ok(())
    }

    #[test]
    fn test_write_data_field_value() -> io::Result<()> {
        use sam::record::data::field::Value;

        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data_field_value(buf, value)?;
            assert_eq!(&buf[..], expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Char('n'), &[b'n'])?;
        t(&mut buf, &Value::Float(8.0), &[0x00, 0x00, 0x00, 0x41])?;

        t(
            &mut buf,
            &Value::String(String::from("ndls")),
            &[b'n', b'd', b'l', b's', 0x00],
        )?;

        t(
            &mut buf,
            &Value::Hex(String::from("cafe")),
            &[b'c', b'a', b'f', b'e', 0x00],
        )?;

        t(
            &mut buf,
            &Value::Int8Array(vec![1, -2]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x01, // values[0] = 1
                0xfe, // values[1] = -2
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt8Array(vec![3, 5]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x03, // values[0] = 3
                0x05, // values[1] = 5
            ],
        )?;

        t(
            &mut buf,
            &Value::Int16Array(vec![8, -13]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x08, 0x00, // values[0] = 8
                0xf3, 0xff, // values[1] = -13
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt16Array(vec![21, 34]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x15, 0x00, // values[0] = 21
                0x22, 0x00, // values[1] = 34
            ],
        )?;

        t(
            &mut buf,
            &Value::Int32Array(vec![55, -89]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x37, 0x00, 0x00, 0x00, // values[0] = 55
                0xa7, 0xff, 0xff, 0xff, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt32Array(vec![144, 223]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x90, 0x00, 0x00, 0x00, // values[0] = 55
                0xdf, 0x00, 0x00, 0x00, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Value::FloatArray(vec![8.0, 13.0]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x00, 0x00, 0x00, 0x41, // values[0] = 8.0
                0x00, 0x00, 0x50, 0x41, // values[1] = 13.0
            ],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_reference_sequence_id() -> Result<(), Box<dyn std::error::Error>> {
        use sam::header::ReferenceSequence;

        let reference_sequences = vec![("sq0", 8), ("sq1", 13)]
            .into_iter()
            .map(|(name, len)| ReferenceSequence::new(name, len).map(|rs| (name.into(), rs)))
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

    #[test]
    fn test_write_position() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let position = sam::record::Position::try_from(8)?;
        write_position(&mut buf, Some(position))?;
        assert_eq!(buf, [0x07, 0x00, 0x00, 0x00]); // pos = 7

        buf.clear();
        write_position(&mut buf, None)?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff]); // pos = -1

        Ok(())
    }

    #[test]
    fn test_region_to_bin() {
        // § 5.3 C source code for computing bin number and overlapping bins (2021-01-07)
        // [-1, 0]
        assert_eq!(region_to_bin(-1, 0), 4680);
        // [8, 13]
        assert_eq!(region_to_bin(7, 13), 4681);
        // [63245986, 63245986]
        assert_eq!(region_to_bin(63245985, 63255986), 8541);
    }
}
