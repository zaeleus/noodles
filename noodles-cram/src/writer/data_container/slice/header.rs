use std::io::{self, Write};

use noodles_core::Position;

use crate::{
    data_container::slice,
    writer::num::{write_itf8, write_ltf8},
};

pub fn write_header<W>(writer: &mut W, header: &slice::Header) -> io::Result<()>
where
    W: Write,
{
    let reference_sequence_id = i32::from(header.reference_sequence_id());
    write_itf8(writer, reference_sequence_id)?;

    write_alignment_start(writer, header.alignment_start())?;

    let alignment_span = i32::try_from(header.alignment_span())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, alignment_span)?;

    let record_count = i32::try_from(header.record_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, record_count)?;

    write_ltf8(writer, header.record_counter())?;

    let block_count = i32::try_from(header.block_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, block_count)?;

    write_block_content_ids(writer, header.block_content_ids())?;

    write_embedded_reference_bases_block_content_id(
        writer,
        header.embedded_reference_bases_block_content_id(),
    )?;

    write_reference_md5(writer, header.reference_md5())?;

    if !header.optional_tags().is_empty() {
        write_optional_tags(writer, header.optional_tags())?;
    }

    Ok(())
}

fn write_alignment_start<W>(writer: &mut W, alignment_start: Option<Position>) -> io::Result<()>
where
    W: Write,
{
    let n = alignment_start.map(usize::from).unwrap_or_default();

    let alignment_start =
        i32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_itf8(writer, alignment_start)
}

fn write_block_content_ids<W>(writer: &mut W, block_content_ids: &[i32]) -> io::Result<()>
where
    W: Write,
{
    let len = block_content_ids.len() as i32;
    write_itf8(writer, len)?;

    for &block_content_id in block_content_ids {
        write_itf8(writer, block_content_id)?;
    }

    Ok(())
}

fn write_embedded_reference_bases_block_content_id<W>(
    writer: &mut W,
    id: Option<i32>,
) -> io::Result<()>
where
    W: Write,
{
    const MISSING: i32 = -1;

    let embedded_reference_bases_block_content_id = id.map(i32::from).unwrap_or(MISSING);
    write_itf8(writer, embedded_reference_bases_block_content_id)
}

fn write_reference_md5<W>(writer: &mut W, reference_md5: &[u8]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(reference_md5)
}

fn write_optional_tags<W>(writer: &mut W, optional_tags: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let len = optional_tags.len() as i32;
    write_itf8(writer, len)?;

    writer.write_all(optional_tags)
}
