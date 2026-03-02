use std::io::{self, Write};

use crate::{
    container::{ReferenceSequenceContext, block, slice},
    file_definition::Version,
    io::writer::num::{write_header_int, write_int, write_long, write_position, write_signed_int},
};

pub fn write_header<W>(writer: &mut W, header: &slice::Header, version: Version) -> io::Result<()>
where
    W: Write,
{
    write_slice_reference_sequence_context(writer, header.reference_sequence_context(), version)?;
    write_record_count(writer, header.record_count(), version)?;
    write_record_counter(writer, header.record_counter(), version)?;
    write_block_count(writer, header.block_count(), version)?;
    write_block_content_ids(writer, header.block_content_ids(), version)?;

    write_embedded_reference_bases_block_content_id(
        writer,
        header.embedded_reference_bases_block_content_id(),
        version,
    )?;

    write_reference_md5(writer, header.reference_md5())?;

    if !header.optional_tags().is_empty() {
        write_optional_tags(writer, header.optional_tags(), version)?;
    }

    Ok(())
}

/// Writes the reference sequence context for a slice header.
///
/// Unlike container headers which use unsigned VLQ (bitcast) for ref_seq_id,
/// slice headers use signed/zigzag encoding for ref_seq_id in CRAM 4.0,
/// matching htslib's `varint_put32s`/`varint_get32s`.
fn write_slice_reference_sequence_context<W>(
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

    // Slice header ref_seq_id uses signed/zigzag encoding in CRAM 4.0,
    // matching htslib's varint_put32s.
    write_signed_int(writer, version, reference_sequence_id)?;
    // alignment_start and alignment_span use 64-bit positions, same as container header.
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

fn write_block_count<W>(writer: &mut W, block_count: usize, version: Version) -> io::Result<()>
where
    W: Write,
{
    let n =
        i32::try_from(block_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, n)
}

fn write_block_content_ids<W>(
    writer: &mut W,
    block_content_ids: &[block::ContentId],
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    let len = i32::try_from(block_content_ids.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, len)?;

    for &block_content_id in block_content_ids {
        write_int(writer, version, block_content_id)?;
    }

    Ok(())
}

fn write_embedded_reference_bases_block_content_id<W>(
    writer: &mut W,
    id: Option<block::ContentId>,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    const MISSING: i32 = -1;

    let embedded_reference_bases_block_content_id = id.unwrap_or(MISSING);
    write_header_int(writer, version, embedded_reference_bases_block_content_id)
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

fn write_optional_tags<W>(writer: &mut W, optional_tags: &[u8], version: Version) -> io::Result<()>
where
    W: Write,
{
    let len = i32::try_from(optional_tags.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, len)?;

    writer.write_all(optional_tags)
}
