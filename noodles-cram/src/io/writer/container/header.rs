use std::io::{self, Write};

use flate2::CrcWriter;

use crate::{
    container::{Header, ReferenceSequenceContext},
    file_definition::Version,
    io::writer::num::{
        write_header_int, write_i32_le, write_int, write_long, write_position, write_u32_le,
    },
};

pub fn write_header<W>(
    writer: &mut W,
    header: &Header,
    len: usize,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    if version.has_crc32() {
        let mut crc_writer = CrcWriter::new(writer);
        write_header_body(&mut crc_writer, header, len, version)?;
        let crc32 = crc_writer.crc().sum();
        write_u32_le(crc_writer.get_mut(), crc32)?;
        Ok(())
    } else {
        write_header_body(writer, header, len, version)
    }
}

fn write_header_body<W>(
    writer: &mut W,
    header: &Header,
    len: usize,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let length = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if version.uses_vlq() {
        write_int(writer, version, length)?;
    } else {
        write_i32_le(writer, length)?;
    }

    write_reference_sequence_context(writer, header.reference_sequence_context(), version)?;
    write_record_count(writer, header.record_count(), version)?;
    write_record_counter(writer, header.record_counter(), version)?;
    write_base_count(writer, header.base_count(), version)?;
    write_block_count(writer, header.block_count(), version)?;
    write_landmarks(writer, header.landmarks(), version)?;

    Ok(())
}

pub(super) fn write_reference_sequence_context<W>(
    writer: &mut W,
    reference_sequence_context: ReferenceSequenceContext,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    const UNMAPPED: i32 = -1;
    const MULTIREF: i32 = -2;

    let (reference_sequence_id, alignment_start, alignment_span) = match reference_sequence_context
    {
        ReferenceSequenceContext::Some(context) => {
            let id = i32::try_from(context.reference_sequence_id())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            let start = i64::try_from(usize::from(context.alignment_start()))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            let span = i64::try_from(context.alignment_span())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            (id, start, span)
        }
        ReferenceSequenceContext::None => (UNMAPPED, 0, 0),
        ReferenceSequenceContext::Many => (MULTIREF, 0, 0),
    };

    // ref_seq_id uses unsigned VLQ (bitcast from i32) in CRAM 4.0, ITF8 in earlier versions.
    write_header_int(writer, version, reference_sequence_id)?;
    // alignment_start and alignment_span use 64-bit positions in CRAM 4.0.
    write_position(writer, version, alignment_start)?;
    write_position(writer, version, alignment_span)?;

    Ok(())
}

fn write_record_count<W>(writer: &mut W, record_count: usize, version: Version) -> io::Result<()>
where
    W: Write,
{
    let n =
        i32::try_from(record_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, n)
}

fn write_record_counter<W>(writer: &mut W, record_counter: u64, version: Version) -> io::Result<()>
where
    W: Write,
{
    let n = i64::try_from(record_counter)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_long(writer, version, n)
}

fn write_base_count<W>(writer: &mut W, base_count: u64, version: Version) -> io::Result<()>
where
    W: Write,
{
    let n =
        i64::try_from(base_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_long(writer, version, n)
}

fn write_block_count<W>(writer: &mut W, block_count: usize, version: Version) -> io::Result<()>
where
    W: Write,
{
    let n =
        i32::try_from(block_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, n)
}

fn write_landmarks<W>(writer: &mut W, landmarks: &[usize], version: Version) -> io::Result<()>
where
    W: Write,
{
    let landmarks_len = i32::try_from(landmarks.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, landmarks_len)?;

    for &pos in landmarks {
        let n = i32::try_from(pos).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_int(writer, version, n)?;
    }

    Ok(())
}
