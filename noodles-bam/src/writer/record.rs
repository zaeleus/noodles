use std::{
    cmp,
    convert::TryFrom,
    ffi::CString,
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};
use noodles_sam::{
    self as sam,
    header::ReferenceSequences,
    record::{Cigar, Data, QualityScores, Sequence},
};

use crate::record::sequence::Base;

// § 4.2 The BAM format (2020-04-30)
//
// ref_id (4) + pos (4) + l_read_name (1) + mapq (1) + bin (2) + n_cigar_op (2) + flag (2) + l_seq
// (4) + next_ref_id (4) + next_pos (4) + tlen (4)
const BLOCK_HEADER_SIZE: usize = 32;

// § 4.2.1 BIN field calculation (2020-04-30)
const UNMAPPED_BIN: u16 = 4680;

// § 4.2.3 SEQ and QUAL encoding (2020-04-30)
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

    let reference_sequence_id =
        find_reference_sequence_id(reference_sequences, record.reference_sequence_name())?;
    let mate_reference_sequence_id =
        find_reference_sequence_id(reference_sequences, record.mate_reference_sequence_name())?;

    let read_name = c_read_name.as_bytes_with_nul();
    let l_read_name = read_name.len() as u8;
    let n_cigar_op = record.cigar().len() as u16;
    let l_seq = record.sequence().len() as i32;
    let data_len = calculate_data_len(record.data()) as i32;

    let block_size = BLOCK_HEADER_SIZE as i32
        + i32::from(l_read_name)
        + (4 * i32::from(n_cigar_op))
        + ((l_seq + 1) / 2)
        + l_seq
        + data_len;

    writer.write_i32::<LittleEndian>(block_size)?;

    let ref_id = reference_sequence_id;
    writer.write_i32::<LittleEndian>(ref_id)?;

    let pos = record
        .position()
        .map(|v| i32::from(v) - 1)
        .unwrap_or(crate::record::UNMAPPED_POSITION);
    writer.write_i32::<LittleEndian>(pos)?;

    writer.write_u8(l_read_name)?;

    let mapq = u8::from(record.mapping_quality());
    writer.write_u8(mapq)?;

    let bin = record
        .position()
        .map(i32::from)
        .map(|start| {
            // 0-based, [start, end)
            let reference_len = record.cigar().reference_len() as i32;
            let end = start + reference_len;
            region_to_bin(start, end) as u16
        })
        .unwrap_or(UNMAPPED_BIN);

    writer.write_u16::<LittleEndian>(bin)?;

    writer.write_u16::<LittleEndian>(n_cigar_op)?;

    let flag = u16::from(record.flags());
    writer.write_u16::<LittleEndian>(flag)?;

    writer.write_i32::<LittleEndian>(l_seq)?;

    let next_ref_id = mate_reference_sequence_id;
    writer.write_i32::<LittleEndian>(next_ref_id)?;

    let next_pos = record
        .mate_position()
        .map(|v| i32::from(v) - 1)
        .unwrap_or(crate::record::UNMAPPED_POSITION);
    writer.write_i32::<LittleEndian>(next_pos)?;

    let tlen = record.template_length();
    writer.write_i32::<LittleEndian>(tlen)?;

    writer.write_all(read_name)?;

    write_cigar(writer, record.cigar())?;

    // § 4.2.3 SEQ and QUAL encoding (2020-04-30)
    let sequence = record.sequence();
    let quality_scores = record.quality_scores();

    write_seq(writer, sequence)?;

    match sequence.len().cmp(&quality_scores.len()) {
        cmp::Ordering::Less => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "quality scores length does not match sequence length",
            ));
        }
        cmp::Ordering::Greater => {
            if quality_scores.is_empty() {
                for _ in 0..sequence.len() {
                    writer.write_u8(NULL_QUALITY_SCORE)?;
                }
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "quality scores length does not match sequence length",
                ));
            }
        }
        cmp::Ordering::Equal => {
            write_qual(writer, quality_scores)?;
        }
    }

    write_data(writer, record.data())?;

    Ok(())
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

        let r = if let Some(c) = chunk.get(1) {
            Base::from(*c)
        } else {
            Base::Eq
        };

        let value = (l as u8) << 4 | (r as u8);

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

fn calculate_data_len(data: &Data) -> usize {
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
            Value::Int32(n) => {
                if *n >= 0 {
                    if *n <= i32::from(u8::MAX) {
                        len += mem::size_of::<u8>();
                    } else if *n <= i32::from(u16::MAX) {
                        len += mem::size_of::<u16>();
                    } else {
                        len += mem::size_of::<u32>();
                    }
                } else if *n >= i32::from(i8::MIN) {
                    len += mem::size_of::<i8>();
                } else if *n >= i32::from(i16::MIN) {
                    len += mem::size_of::<i16>();
                } else {
                    len += mem::size_of::<i32>();
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

    len
}

fn write_data<W>(writer: &mut W, data: &Data) -> io::Result<()>
where
    W: Write,
{
    use noodles_sam::record::data::field::Value;

    for field in data.values() {
        writer.write_all(field.tag().as_ref().as_bytes())?;

        let value = field.value();

        if let Value::Int32(n) = value {
            write_data_i32_value(writer, *n)?;
            continue;
        }

        writer.write_u8(char::from(value.ty()) as u8)?;

        if let Some(subtype) = value.subtype() {
            writer.write_u8(char::from(subtype) as u8)?;
        }

        match value {
            Value::Char(c) => {
                writer.write_u8(*c as u8)?;
            }
            Value::Int32(_) => unreachable!(),
            Value::Float(n) => {
                writer.write_f32::<LittleEndian>(*n)?;
            }
            Value::String(s) | Value::Hex(s) => {
                let c_str = CString::new(s.as_bytes())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                writer.write_all(c_str.as_bytes_with_nul())?;
            }
            Value::Int8Array(values) => {
                writer.write_u32::<LittleEndian>(values.len() as u32)?;

                for &n in values {
                    writer.write_i8(n)?;
                }
            }
            Value::UInt8Array(values) => {
                writer.write_u32::<LittleEndian>(values.len() as u32)?;

                for &n in values {
                    writer.write_u8(n)?;
                }
            }
            Value::Int16Array(values) => {
                writer.write_u32::<LittleEndian>(values.len() as u32)?;

                for &n in values {
                    writer.write_i16::<LittleEndian>(n)?;
                }
            }
            Value::UInt16Array(values) => {
                writer.write_u32::<LittleEndian>(values.len() as u32)?;

                for &n in values {
                    writer.write_u16::<LittleEndian>(n)?;
                }
            }
            Value::Int32Array(values) => {
                writer.write_u32::<LittleEndian>(values.len() as u32)?;

                for &n in values {
                    writer.write_i32::<LittleEndian>(n)?;
                }
            }
            Value::UInt32Array(values) => {
                writer.write_u32::<LittleEndian>(values.len() as u32)?;

                for &n in values {
                    writer.write_u32::<LittleEndian>(n)?;
                }
            }
            Value::FloatArray(values) => {
                writer.write_u32::<LittleEndian>(values.len() as u32)?;

                for &n in values {
                    writer.write_f32::<LittleEndian>(n)?;
                }
            }
        }
    }

    Ok(())
}

fn write_data_i32_value<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    use crate::record::data::field::value::Type;

    if n >= 0 {
        if n <= i32::from(u8::MAX) {
            writer.write_u8(u8::from(Type::UInt8))?;
            writer.write_u8(n as u8)
        } else if n <= i32::from(u16::MAX) {
            writer.write_u8(u8::from(Type::UInt16))?;
            writer.write_u16::<LittleEndian>(n as u16)
        } else {
            writer.write_u8(u8::from(Type::UInt32))?;
            writer.write_u32::<LittleEndian>(n as u32)
        }
    } else if n >= i32::from(i8::MIN) {
        writer.write_u8(u8::from(Type::Int8))?;
        writer.write_i8(n as i8)
    } else if n >= i32::from(i16::MIN) {
        writer.write_u8(u8::from(Type::Int16))?;
        writer.write_i16::<LittleEndian>(n as i16)
    } else {
        writer.write_u8(u8::from(Type::Int32))?;
        writer.write_i32::<LittleEndian>(n)
    }
}

fn find_reference_sequence_id(
    reference_sequences: &ReferenceSequences,
    reference_sequence_name: Option<&sam::record::ReferenceSequenceName>,
) -> io::Result<i32> {
    use crate::record::reference_sequence_id;

    match reference_sequence_name {
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
            }),
        None => Ok(reference_sequence_id::UNMAPPED),
    }
}

// § 5.3 C source code for computing bin number and overlapping bins (2020-04-30)
// 0-based, [start, end)
#[allow(clippy::eq_op)]
fn region_to_bin(start: i32, mut end: i32) -> i32 {
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
    fn test_write_data_i32_value() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: i32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data_i32_value(buf, n)?;
            assert_eq!(&buf[..], expected, "n = {}", n);
            Ok(())
        }

        let mut buf = Vec::new();

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

        Ok(())
    }

    #[test]
    fn test_find_reference_sequence_id() -> Result<(), Box<dyn std::error::Error>> {
        use sam::header::ReferenceSequence;

        let reference_sequences = vec![("sq0", 8), ("sq1", 13), ("sq2", 21)]
            .into_iter()
            .map(|(name, len)| ReferenceSequence::new(name, len).map(|rs| (name.into(), rs)))
            .collect::<Result<_, _>>()?;

        assert!(matches!(
            find_reference_sequence_id(&reference_sequences, None),
            Ok(-1)
        ));

        let name = "sq1".parse().map(Some)?;
        assert!(matches!(
            find_reference_sequence_id(&reference_sequences, name.as_ref()),
            Ok(1)
        ));

        let name = "sq3".parse().map(Some)?;
        assert!(matches!(
            find_reference_sequence_id(&reference_sequences, name.as_ref()),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

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
