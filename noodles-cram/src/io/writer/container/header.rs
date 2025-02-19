use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;

use crate::{
    container::{Header, ReferenceSequenceContext},
    io::writer::num::{write_itf8, write_ltf8},
};

pub fn write_header<W>(writer: &mut W, header: &Header, len: usize) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    let length = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    crc_writer.write_i32::<LittleEndian>(length)?;

    write_reference_sequence_context(&mut crc_writer, header.reference_sequence_context())?;

    let number_of_records = i32::try_from(header.record_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, number_of_records)?;

    let record_counter = i64::try_from(header.record_counter())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_ltf8(&mut crc_writer, record_counter)?;

    let bases = i64::try_from(header.base_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_ltf8(&mut crc_writer, bases)?;

    let number_of_blocks = i32::try_from(header.block_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, number_of_blocks)?;

    write_landmarks(&mut crc_writer, header.landmarks())?;

    let crc32 = crc_writer.crc().sum();
    let writer = crc_writer.into_inner();
    writer.write_u32::<LittleEndian>(crc32)?;

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
