use std::io::{self, Write};

use flate2::CrcWriter;

use crate::{
    container::{Header, ReferenceSequenceContext},
    io::writer::num::{write_i32_le, write_itf8, write_ltf8, write_u32_le},
};

pub fn write_header<W>(writer: &mut W, header: &Header, len: usize) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);
    write_header_inner(&mut crc_writer, header, len)
}

fn write_header_inner<W>(writer: &mut CrcWriter<W>, header: &Header, len: usize) -> io::Result<()>
where
    W: Write,
{
    let length = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, length)?;

    write_reference_sequence_context(writer, header.reference_sequence_context())?;
    write_record_count(writer, header.record_count())?;
    write_record_counter(writer, header.record_counter())?;
    write_base_count(writer, header.base_count())?;
    write_block_count(writer, header.block_count())?;
    write_landmarks(writer, header.landmarks())?;

    let crc32 = writer.crc().sum();
    write_u32_le(writer.get_mut(), crc32)?;

    Ok(())
}

pub(super) fn write_reference_sequence_context<W>(
    writer: &mut W,
    reference_sequence_context: ReferenceSequenceContext,
) -> io::Result<()>
where
    W: Write,
{
    const MISSING: i32 = 0;
    const UNMAPPED: i32 = -1;
    const MULTIREF: i32 = -2;

    let (reference_sequence_id, alignment_start, alignment_span) = match reference_sequence_context
    {
        ReferenceSequenceContext::Some(context) => {
            let id = i32::try_from(context.reference_sequence_id())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            let start = i32::try_from(usize::from(context.alignment_start()))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            let span = i32::try_from(context.alignment_span())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            (id, start, span)
        }
        ReferenceSequenceContext::None => (UNMAPPED, MISSING, MISSING),
        ReferenceSequenceContext::Many => (MULTIREF, MISSING, MISSING),
    };

    write_itf8(writer, reference_sequence_id)?;
    write_itf8(writer, alignment_start)?;
    write_itf8(writer, alignment_span)?;

    Ok(())
}

fn write_record_count<W>(writer: &mut W, record_count: usize) -> io::Result<()>
where
    W: Write,
{
    let n =
        i32::try_from(record_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, n)
}

fn write_record_counter<W>(writer: &mut W, record_counter: u64) -> io::Result<()>
where
    W: Write,
{
    let n = i64::try_from(record_counter)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_ltf8(writer, n)
}

fn write_base_count<W>(writer: &mut W, base_count: u64) -> io::Result<()>
where
    W: Write,
{
    let n =
        i64::try_from(base_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_ltf8(writer, n)
}

fn write_block_count<W>(writer: &mut W, block_count: usize) -> io::Result<()>
where
    W: Write,
{
    let n =
        i32::try_from(block_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, n)
}

fn write_landmarks<W>(writer: &mut W, landmarks: &[usize]) -> io::Result<()>
where
    W: Write,
{
    let landmarks_len = i32::try_from(landmarks.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, landmarks_len)?;

    for &pos in landmarks {
        let n = i32::try_from(pos).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_itf8(writer, n)?;
    }

    Ok(())
}
