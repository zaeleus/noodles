use std::{
    convert::TryFrom,
    io::{self, Write},
};

use crate::{
    container::slice,
    num::Itf8,
    writer::num::{write_itf8, write_ltf8},
};

pub fn write_header<W>(writer: &mut W, header: &slice::Header) -> io::Result<()>
where
    W: Write,
{
    let reference_sequence_id = i32::from(header.reference_sequence_id());
    write_itf8(writer, reference_sequence_id)?;

    write_itf8(writer, header.alignment_start())?;
    write_itf8(writer, header.alignment_span())?;
    write_itf8(writer, header.record_count())?;
    write_ltf8(writer, header.record_counter())?;

    let block_count = Itf8::try_from(header.block_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, block_count)?;

    write_block_content_ids(writer, header.block_content_ids())?;

    let embedded_reference_bases_block_content_id =
        i32::from(header.embedded_reference_bases_block_content_id());
    write_itf8(writer, embedded_reference_bases_block_content_id)?;

    write_reference_md5(writer, header.reference_md5())?;

    if !header.optional_tags().is_empty() {
        write_optional_tags(writer, header.optional_tags())?;
    }

    Ok(())
}

fn write_block_content_ids<W>(writer: &mut W, block_content_ids: &[Itf8]) -> io::Result<()>
where
    W: Write,
{
    let len = block_content_ids.len() as Itf8;
    write_itf8(writer, len)?;

    for &block_content_id in block_content_ids {
        write_itf8(writer, block_content_id)?;
    }

    Ok(())
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
    let len = optional_tags.len() as Itf8;
    write_itf8(writer, len)?;

    writer.write_all(optional_tags)
}
