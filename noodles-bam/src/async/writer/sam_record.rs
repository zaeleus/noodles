use std::{ffi::CString, num};

use noodles_sam::{self as sam, header::ReferenceSequences};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

// § 1.4 "The alignment section: mandatory fields" (2021-06-03): "A `QNAME` '*' indicates the
// information is unavailable."
const MISSING_READ_NAME: &str = "*";

pub async fn write_sam_record<W>(
    writer: &mut W,
    reference_sequences: &ReferenceSequences,
    record: &sam::Record,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::writer::sam_record::calculate_data_len;

    let name = record
        .read_name()
        .map(|name| name.as_str())
        .unwrap_or(MISSING_READ_NAME);
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

    let block_size = 32
        + u32::from(l_read_name)
        + (4 * u32::from(n_cigar_op))
        + ((l_seq + 1) / 2)
        + l_seq
        + data_len;
    writer.write_u32_le(block_size).await?;

    // ref_id
    write_reference_sequence_id(
        writer,
        reference_sequences,
        record.reference_sequence_name(),
    )
    .await?;

    // pos
    write_position(writer, record.position()).await?;

    // l_read_name
    writer.write_u8(l_read_name).await?;

    // mapq
    write_mapping_quality(writer, record.mapping_quality()).await?;

    // bin
    write_bin(writer, record).await?;

    // n_cigar_op
    writer.write_u16_le(n_cigar_op).await?;

    // flag
    write_flags(writer, record.flags()).await?;

    // l_seq
    writer.write_u32_le(l_seq).await?;

    // next_ref_id
    write_reference_sequence_id(
        writer,
        reference_sequences,
        record.mate_reference_sequence_name(),
    )
    .await?;

    // next_pos
    write_position(writer, record.mate_position()).await?;

    // tlen
    write_template_length(writer, record.template_length()).await?;

    // read_name
    writer.write_all(read_name).await?;

    // cigar
    write_cigar(writer, record.cigar()).await?;

    let sequence = record.sequence();

    if !sequence.is_empty() {
        // seq
        write_sequence(writer, record.sequence()).await?;

        let quality_scores = record.quality_scores();

        // qual
        if sequence.len() == quality_scores.len() {
            write_quality_scores(writer, record.quality_scores()).await?;
        } else if quality_scores.is_empty() {
            write_missing_filled_quality_scores(writer, sequence.len()).await?;
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

    write_data(writer, record.data()).await?;

    Ok(())
}

async fn write_reference_sequence_id<W>(
    writer: &mut W,
    reference_sequences: &ReferenceSequences,
    reference_sequence_name: Option<&sam::record::ReferenceSequenceName>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
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

    writer.write_i32_le(id).await
}

async fn write_position<W>(
    writer: &mut W,
    position: Option<sam::record::Position>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::record::UNMAPPED_POSITION;

    let pos = position
        .map(|p| i32::from(p) - 1)
        .unwrap_or(UNMAPPED_POSITION);

    writer.write_i32_le(pos).await
}

async fn write_mapping_quality<W>(
    writer: &mut W,
    mapping_quality: sam::record::MappingQuality,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let mapq = u8::from(mapping_quality);
    writer.write_u8(mapq).await
}

async fn write_bin<W>(writer: &mut W, record: &sam::Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::writer::sam_record::region_to_bin;

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

    writer.write_u16_le(bin).await
}

async fn write_flags<W>(writer: &mut W, flags: sam::record::Flags) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let flag = u16::from(flags);
    writer.write_u16_le(flag).await
}

async fn write_template_length<W>(writer: &mut W, tlen: i32) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_i32_le(tlen).await
}

async fn write_cigar<W>(writer: &mut W, cigar: &sam::record::Cigar) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    for op in cigar.iter() {
        let len = op.len();
        let kind = op.kind() as u32;
        let value = len << 4 | kind;
        writer.write_u32_le(value).await?;
    }

    Ok(())
}

async fn write_sequence<W>(writer: &mut W, sequence: &sam::record::Sequence) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::record::sequence::Base;

    for chunk in sequence.chunks(2) {
        let l = Base::from(chunk[0]);

        // § 4.2.3 "SEQ and QUAL encoding" (2021-06-03): "When `l_seq` is odd the bottom 4 bits of
        // the last byte are undefined, but we recommend writing these are zero."
        let r = chunk.get(1).copied().map(Base::from).unwrap_or(Base::Eq);

        let value = u8::from(l) << 4 | u8::from(r);
        writer.write_u8(value).await?;
    }

    Ok(())
}

async fn write_quality_scores<W>(
    writer: &mut W,
    quality_scores: &sam::record::QualityScores,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    for &score in quality_scores.iter() {
        let value = u8::from(score);
        writer.write_u8(value).await?;
    }

    Ok(())
}

async fn write_missing_filled_quality_scores<W>(
    writer: &mut W,
    sequence_len: usize,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    // § 4.2.3 "SEQ and QUAL encoding" (2021-06-03): "When base qualities are omitted but the
    // sequence is not, `qual` is filled with `OxFF` bytes (to length `l_seq`)."
    const MISSING_QUALITY_SCORE: u8 = 0xff;

    for _ in 0..sequence_len {
        writer.write_u8(MISSING_QUALITY_SCORE).await?;
    }

    Ok(())
}

async fn write_data<W>(writer: &mut W, data: &sam::record::Data) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    for field in data.values() {
        write_data_field(writer, field).await?;
    }

    Ok(())
}

async fn write_data_field<W>(writer: &mut W, field: &sam::record::data::Field) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use sam::record::data::field::Value;

    write_data_field_tag(writer, field.tag()).await?;

    let value = field.value();

    if let Value::Int(n) = value {
        write_data_field_int_value(writer, *n).await?;
    } else {
        write_data_field_value_type(writer, value).await?;
        write_data_field_value(writer, value).await?;
    }

    Ok(())
}

async fn write_data_field_tag<W>(
    writer: &mut W,
    tag: sam::record::data::field::Tag,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_all(tag.as_ref()).await
}

async fn write_data_field_int_value<W>(writer: &mut W, n: i64) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::record::data::field::value::Type;

    if n >= 0 {
        if n <= i64::from(u8::MAX) {
            writer.write_u8(u8::from(Type::UInt8)).await?;
            return writer.write_u8(n as u8).await;
        } else if n <= i64::from(u16::MAX) {
            writer.write_u8(u8::from(Type::UInt16)).await?;
            return writer.write_u16_le(n as u16).await;
        } else if n <= i64::from(u32::MAX) {
            writer.write_u8(u8::from(Type::UInt32)).await?;
            return writer.write_u32_le(n as u32).await;
        }
    } else if n >= i64::from(i8::MIN) {
        writer.write_u8(u8::from(Type::Int8)).await?;
        return writer.write_i8(n as i8).await;
    } else if n >= i64::from(i16::MIN) {
        writer.write_u8(u8::from(Type::Int16)).await?;
        return writer.write_i16_le(n as i16).await;
    } else if n >= i64::from(i32::MIN) {
        writer.write_u8(u8::from(Type::Int32)).await?;
        return writer.write_i32_le(n as i32).await;
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("invalid data field integer value: {}", n),
    ))
}

async fn write_data_field_value_type<W>(
    writer: &mut W,
    value: &sam::record::data::field::Value,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let val_type = char::from(value.ty()) as u8;
    writer.write_u8(val_type).await?;

    if let Some(subtype) = value.subtype() {
        let val_subtype = char::from(subtype) as u8;
        writer.write_u8(val_subtype).await?;
    }

    Ok(())
}

async fn write_data_field_value<W>(
    writer: &mut W,
    value: &sam::record::data::field::Value,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use sam::record::data::field::Value;

    fn invalid_array_len(e: num::TryFromIntError) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidInput, e)
    }

    match value {
        Value::Char(c) => writer.write_u8(*c as u8).await?,
        Value::Int(_) => {
            // Integers are handled by `write_data_field_int_value`.
            unreachable!();
        }
        Value::Float(n) => writer.write_f32_le(*n).await?,
        Value::String(s) | Value::Hex(s) => {
            let c_str = CString::new(s.as_bytes())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            writer.write_all(c_str.as_bytes_with_nul()).await?;
        }
        Value::Int8Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32_le(len).await?;

            for &n in values {
                writer.write_i8(n).await?;
            }
        }
        Value::UInt8Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32_le(len).await?;

            for &n in values {
                writer.write_u8(n).await?;
            }
        }
        Value::Int16Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32_le(len).await?;

            for &n in values {
                writer.write_i16_le(n).await?;
            }
        }
        Value::UInt16Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32_le(len).await?;

            for &n in values {
                writer.write_u16_le(n).await?;
            }
        }
        Value::Int32Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32_le(len).await?;

            for &n in values {
                writer.write_i32_le(n).await?;
            }
        }
        Value::UInt32Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32_le(len).await?;

            for &n in values {
                writer.write_u32_le(n).await?;
            }
        }
        Value::FloatArray(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            writer.write_u32_le(len).await?;

            for &n in values {
                writer.write_f32_le(n).await?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_reference_sequence_id() -> Result<(), Box<dyn std::error::Error>> {
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
        )
        .await?;
        assert_eq!(buf, [0x00, 0x00, 0x00, 0x00]);

        buf.clear();
        write_reference_sequence_id(&mut buf, &reference_sequences, None).await?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff]);

        buf.clear();
        let reference_sequence_name = "sq2".parse()?;
        assert!(matches!(
            write_reference_sequence_id(&mut buf, &reference_sequences, Some(&reference_sequence_name)).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[tokio::test]
    async fn test_write_position() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let position = sam::record::Position::try_from(8)?;
        write_position(&mut buf, Some(position)).await?;
        assert_eq!(buf, [0x07, 0x00, 0x00, 0x00]); // pos = 7

        buf.clear();
        write_position(&mut buf, None).await?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff]); // pos = -1

        Ok(())
    }

    #[tokio::test]
    async fn test_write_missing_filled_quality_scores() -> io::Result<()> {
        let mut buf = Vec::new();
        write_missing_filled_quality_scores(&mut buf, 4).await?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff]);
        Ok(())
    }

    #[tokio::test]
    async fn test_write_data_field_int_value() -> io::Result<()> {
        async fn t(buf: &mut Vec<u8>, n: i64, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data_field_int_value(buf, n).await?;
            assert_eq!(&buf[..], expected);
            Ok(())
        }

        let mut buf = Vec::new();

        // i32::MIN - 1
        buf.clear();
        assert!(matches!(
            write_data_field_int_value(&mut buf, -2147483649).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));
        // i32::MIN
        t(&mut buf, -2147483648, &[b'i', 0x00, 0x00, 0x00, 0x80]).await?;
        // i32::MIN + 1
        t(&mut buf, -2147483647, &[b'i', 0x01, 0x00, 0x00, 0x80]).await?;

        // i16::MIN - 1
        t(&mut buf, -32769, &[b'i', 0xff, 0x7f, 0xff, 0xff]).await?;
        // i16::MIN
        t(&mut buf, -32768, &[b's', 0x00, 0x80]).await?;
        // i16::MIN + 1
        t(&mut buf, -32767, &[b's', 0x01, 0x80]).await?;

        // i8::MIN - 1
        t(&mut buf, -129, &[b's', 0x7f, 0xff]).await?;
        // i8::MIN
        t(&mut buf, -128, &[b'c', 0x80]).await?;
        // i8::MIN + 1
        t(&mut buf, -127, &[b'c', 0x81]).await?;

        // -1
        t(&mut buf, -1, &[b'c', 0xff]).await?;
        // 0
        t(&mut buf, 0, &[b'C', 0x00]).await?;
        // 1
        t(&mut buf, 1, &[b'C', 0x01]).await?;

        // i8::MAX - 1
        t(&mut buf, 126, &[b'C', 0x7e]).await?;
        // i8::MAX
        t(&mut buf, 127, &[b'C', 0x7f]).await?;
        // i8::MAX + 1
        t(&mut buf, 128, &[b'C', 0x80]).await?;

        // u8::MAX - 1
        t(&mut buf, 254, &[b'C', 0xfe]).await?;
        // u8::MAX
        t(&mut buf, 255, &[b'C', 0xff]).await?;
        // u8::MAX + 1
        t(&mut buf, 256, &[b'S', 0x00, 0x01]).await?;

        // i16::MAX - 1
        t(&mut buf, 32766, &[b'S', 0xfe, 0x7f]).await?;
        // i16::MAX
        t(&mut buf, 32767, &[b'S', 0xff, 0x7f]).await?;
        // i16::MAX + 1
        t(&mut buf, 32768, &[b'S', 0x00, 0x80]).await?;

        // u16::MAX - 1
        t(&mut buf, 65534, &[b'S', 0xfe, 0xff]).await?;
        // u16::MAX
        t(&mut buf, 65535, &[b'S', 0xff, 0xff]).await?;
        // u16::MAX + 1
        t(&mut buf, 65536, &[b'I', 0x00, 0x00, 0x01, 0x00]).await?;

        // i32::MAX - 1
        t(&mut buf, 2147483646, &[b'I', 0xfe, 0xff, 0xff, 0x7f]).await?;
        // i32::MAX
        t(&mut buf, 2147483647, &[b'I', 0xff, 0xff, 0xff, 0x7f]).await?;
        // i32::MAX + 1
        t(&mut buf, 2147483648, &[b'I', 0x00, 0x00, 0x00, 0x80]).await?;

        // u32::MAX - 1
        t(&mut buf, 4294967295, &[b'I', 0xff, 0xff, 0xff, 0xff]).await?;
        // u32::MAX
        buf.clear();
        assert!(matches!(
            write_data_field_int_value(&mut buf, 4294967296).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[tokio::test]
    async fn test_write_data_field_value_type() -> io::Result<()> {
        use sam::record::data::field::Value;

        async fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data_field_value_type(buf, value).await?;
            assert_eq!(&buf[..], expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Char('n'), &[b'A']).await?;
        t(&mut buf, &Value::Int(13), &[b'i']).await?;
        t(&mut buf, &Value::Float(8.0), &[b'f']).await?;
        t(&mut buf, &Value::String(String::from("ndls")), &[b'Z']).await?;
        t(&mut buf, &Value::Hex(String::from("cafe")), &[b'H']).await?;
        t(&mut buf, &Value::Int8Array(vec![1, -2]), &[b'B', b'c']).await?;
        t(&mut buf, &Value::UInt8Array(vec![3, 5]), &[b'B', b'C']).await?;
        t(&mut buf, &Value::Int16Array(vec![8, -13]), &[b'B', b's']).await?;
        t(&mut buf, &Value::UInt16Array(vec![21, 34]), &[b'B', b'S']).await?;
        t(&mut buf, &Value::Int32Array(vec![55, -89]), &[b'B', b'i']).await?;
        t(&mut buf, &Value::UInt32Array(vec![144, 223]), &[b'B', b'I']).await?;
        t(&mut buf, &Value::FloatArray(vec![8.0, 13.0]), &[b'B', b'f']).await?;

        Ok(())
    }

    #[tokio::test]
    async fn test_write_data_field_value() -> io::Result<()> {
        use sam::record::data::field::Value;

        async fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data_field_value(buf, value).await?;
            assert_eq!(&buf[..], expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Value::Char('n'), &[b'n']).await?;
        t(&mut buf, &Value::Float(8.0), &[0x00, 0x00, 0x00, 0x41]).await?;

        t(
            &mut buf,
            &Value::String(String::from("ndls")),
            &[b'n', b'd', b'l', b's', 0x00],
        )
        .await?;

        t(
            &mut buf,
            &Value::Hex(String::from("cafe")),
            &[b'c', b'a', b'f', b'e', 0x00],
        )
        .await?;

        t(
            &mut buf,
            &Value::Int8Array(vec![1, -2]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x01, // values[0] = 1
                0xfe, // values[1] = -2
            ],
        )
        .await?;

        t(
            &mut buf,
            &Value::UInt8Array(vec![3, 5]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x03, // values[0] = 3
                0x05, // values[1] = 5
            ],
        )
        .await?;

        t(
            &mut buf,
            &Value::Int16Array(vec![8, -13]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x08, 0x00, // values[0] = 8
                0xf3, 0xff, // values[1] = -13
            ],
        )
        .await?;

        t(
            &mut buf,
            &Value::UInt16Array(vec![21, 34]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x15, 0x00, // values[0] = 21
                0x22, 0x00, // values[1] = 34
            ],
        )
        .await?;

        t(
            &mut buf,
            &Value::Int32Array(vec![55, -89]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x37, 0x00, 0x00, 0x00, // values[0] = 55
                0xa7, 0xff, 0xff, 0xff, // values[1] = -89
            ],
        )
        .await?;

        t(
            &mut buf,
            &Value::UInt32Array(vec![144, 223]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x90, 0x00, 0x00, 0x00, // values[0] = 55
                0xdf, 0x00, 0x00, 0x00, // values[1] = -89
            ],
        )
        .await?;

        t(
            &mut buf,
            &Value::FloatArray(vec![8.0, 13.0]),
            &[
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x00, 0x00, 0x00, 0x41, // values[0] = 8.0
                0x00, 0x00, 0x50, 0x41, // values[1] = 13.0
            ],
        )
        .await?;

        Ok(())
    }
}
