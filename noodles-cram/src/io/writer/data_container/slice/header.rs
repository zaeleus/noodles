use std::io::{self, Write};

use crate::{
    container::block,
    data_container::{slice, ReferenceSequenceContext},
    io::writer::num::{write_itf8, write_ltf8},
};

pub fn write_header<W>(writer: &mut W, header: &slice::Header) -> io::Result<()>
where
    W: Write,
{
    write_reference_sequence_context(writer, header.reference_sequence_context())?;

    let record_count = i32::try_from(header.record_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, record_count)?;

    let record_counter = i64::try_from(header.record_counter())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_ltf8(writer, record_counter)?;

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

fn write_reference_sequence_context<W>(
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
    let len = i32::try_from(optional_tags.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, len)?;

    writer.write_all(optional_tags)
}
