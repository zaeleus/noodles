use std::{ffi::CString, num};

use noodles_sam::{self as sam, header::ReferenceSequences, AlignmentRecord};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use super::record::{
    write_bin, write_cigar, write_flags, write_mapping_quality, write_position,
    write_quality_scores, write_sequence, write_template_length,
};

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
        .map(|name| name.as_ref())
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

    let data_len = calculate_data_len(record.data()).and_then(|n| {
        u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
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
    write_position(writer, record.alignment_start()).await?;

    // l_read_name
    writer.write_u8(l_read_name).await?;

    // mapq
    write_mapping_quality(writer, record.mapping_quality()).await?;

    // bin
    write_bin(writer, record.alignment_start(), record.alignment_end()).await?;

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
    write_position(writer, record.mate_alignment_start()).await?;

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
    write_data_field_tag(writer, field.tag()).await?;
    write_data_field_value_type(writer, field.value()).await?;
    write_data_field_value(writer, field.value()).await?;
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

async fn write_data_field_value_type<W>(
    writer: &mut W,
    value: &sam::record::data::field::Value,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let val_type = u8::from(value.ty());
    writer.write_u8(val_type).await?;

    if let Some(subtype) = value.subtype() {
        let val_subtype = u8::from(subtype);
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
        Value::Int8(n) => writer.write_i8(*n).await?,
        Value::UInt8(n) => writer.write_u8(*n).await?,
        Value::Int16(n) => writer.write_i16_le(*n).await?,
        Value::UInt16(n) => writer.write_u16_le(*n).await?,
        Value::Int32(n) => writer.write_i32_le(*n).await?,
        Value::UInt32(n) => writer.write_u32_le(*n).await?,
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
    async fn test_write_missing_filled_quality_scores() -> io::Result<()> {
        let mut buf = Vec::new();
        write_missing_filled_quality_scores(&mut buf, 4).await?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff]);
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
        t(&mut buf, &Value::Int8(0), &[b'c']).await?;
        t(&mut buf, &Value::UInt8(0), &[b'C']).await?;
        t(&mut buf, &Value::Int16(0), &[b's']).await?;
        t(&mut buf, &Value::UInt16(0), &[b'S']).await?;
        t(&mut buf, &Value::Int32(0), &[b'i']).await?;
        t(&mut buf, &Value::UInt32(0), &[b'I']).await?;
        t(&mut buf, &Value::Float(0.0), &[b'f']).await?;
        t(&mut buf, &Value::String(String::from("ndls")), &[b'Z']).await?;
        t(&mut buf, &Value::Hex(String::from("CAFE")), &[b'H']).await?;
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
        t(&mut buf, &Value::Int8(1), &[0x01]).await?;
        t(&mut buf, &Value::UInt8(2), &[0x02]).await?;
        t(&mut buf, &Value::Int16(3), &[0x03, 0x00]).await?;
        t(&mut buf, &Value::UInt16(5), &[0x05, 0x00]).await?;
        t(&mut buf, &Value::Int32(8), &[0x08, 0x00, 0x00, 0x00]).await?;
        t(&mut buf, &Value::UInt32(13), &[0x0d, 0x00, 0x00, 0x00]).await?;
        t(&mut buf, &Value::Float(8.0), &[0x00, 0x00, 0x00, 0x41]).await?;

        t(
            &mut buf,
            &Value::String(String::from("ndls")),
            &[b'n', b'd', b'l', b's', 0x00],
        )
        .await?;

        t(
            &mut buf,
            &Value::Hex(String::from("CAFE")),
            &[b'C', b'A', b'F', b'E', 0x00],
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
