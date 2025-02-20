use std::io::{self, Write};

use crate::{
    container::{block, slice},
    io::writer::{
        container::header::write_reference_sequence_context,
        num::{write_itf8, write_ltf8},
    },
};

pub fn write_header<W>(writer: &mut W, header: &slice::Header) -> io::Result<()>
where
    W: Write,
{
    write_reference_sequence_context(writer, header.reference_sequence_context())?;
    write_record_count(writer, header.record_count())?;
    write_record_counter(writer, header.record_counter())?;
    write_block_count(writer, header.block_count())?;
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

fn write_block_count<W>(writer: &mut W, block_count: usize) -> io::Result<()>
where
    W: Write,
{
    let n =
        i32::try_from(block_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, n)
}

fn write_block_content_ids<W>(
    writer: &mut W,
    block_content_ids: &[block::ContentId],
) -> io::Result<()>
where
    W: Write,
{
    let len = i32::try_from(block_content_ids.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, len)?;

    for &block_content_id in block_content_ids {
        write_itf8(writer, block_content_id)?;
    }

    Ok(())
}

fn write_embedded_reference_bases_block_content_id<W>(
    writer: &mut W,
    id: Option<block::ContentId>,
) -> io::Result<()>
where
    W: Write,
{
    const MISSING: i32 = -1;

    let embedded_reference_bases_block_content_id = id.unwrap_or(MISSING);
    write_itf8(writer, embedded_reference_bases_block_content_id)
}

fn write_reference_md5<W>(writer: &mut W, md5: Option<&[u8; 16]>) -> io::Result<()>
where
    W: Write,
{
    const MISSING: [u8; 16] = [
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00,
    ];

    let buf = md5.unwrap_or(&MISSING);
    writer.write_all(buf)
}

fn write_optional_tags<W>(writer: &mut W, optional_tags: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let len = i32::try_from(optional_tags.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, len)?;

    writer.write_all(optional_tags)
}
