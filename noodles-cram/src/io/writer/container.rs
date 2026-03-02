mod block;
pub(crate) mod compression_header;
mod header;
pub(crate) mod slice;

use std::{
    cmp,
    collections::HashSet,
    io::{self, Write},
};

use noodles_fasta as fasta;
use noodles_sam as sam;

pub use self::{
    block::{Block, write_block},
    header::write_header,
};
use self::{
    compression_header::{build_compression_header, write_compression_header},
    slice::{Slice, build_slice},
};
use super::{Options, Record};
use crate::{
    container::{Header, ReferenceSequenceContext, block::ContentType},
    file_definition::Version,
};

pub fn write_container<W>(
    writer: &mut W,
    reference_sequence_repository: &fasta::Repository,
    options: &Options,
    header: &sam::Header,
    record_counter: u64,
    records: &mut [Record],
) -> io::Result<()>
where
    W: Write,
{
    if records.is_empty() {
        return Ok(());
    }

    let version = options.version;

    let (header, container_size, blocks) = build_container(
        reference_sequence_repository,
        options,
        header,
        record_counter,
        records,
    )?;

    write_header(writer, &header, container_size, version)?;

    for block in blocks {
        write_block(writer, &block, version)?;
    }

    Ok(())
}

fn build_container(
    reference_sequence_repository: &fasta::Repository,
    options: &Options,
    header: &sam::Header,
    record_counter: u64,
    records: &mut [Record],
) -> io::Result<(Header, usize, Vec<Block>)> {
    let version = options.version;

    let compression_header = build_compression_header(options, records);

    let records_per_slice = options.records_per_slice;

    #[cfg(feature = "parallel")]
    let slices = {
        use rayon::prelude::*;

        records
            .par_chunks_mut(records_per_slice)
            .enumerate()
            .map(|(i, chunk)| {
                let slice_record_counter = record_counter + (i * records_per_slice) as u64;

                build_slice(
                    reference_sequence_repository,
                    options,
                    header,
                    slice_record_counter,
                    &compression_header,
                    chunk,
                )
            })
            .collect::<io::Result<Vec<_>>>()?
    };

    #[cfg(not(feature = "parallel"))]
    let slices = {
        let mut slices = Vec::new();
        let mut slice_record_counter = record_counter;

        for chunk in records.chunks_mut(records_per_slice) {
            let slice = build_slice(
                reference_sequence_repository,
                options,
                header,
                slice_record_counter,
                &compression_header,
                chunk,
            )?;

            slices.push(slice);

            let record_count = u64::try_from(chunk.len())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            slice_record_counter += record_count;
        }

        slices
    };

    let reference_sequence_context = get_container_reference_sequence_context(&slices)?;
    let record_count = records.len();
    let base_count = calculate_base_count(records)?;

    // Collect used block content IDs from all slices to prune unused data
    // series from the compression header. This prevents declaring encodings
    // for data series that have no corresponding blocks.
    let used_content_ids: HashSet<_> = slices
        .iter()
        .flat_map(|s| s.external_data_blocks.iter().map(|b| b.content_id))
        .collect();

    let mut pruned_compression_header = compression_header;
    pruned_compression_header
        .data_series_encodings
        .retain_used_content_ids(&used_content_ids);

    let mut buf = Vec::new();

    write_compression_header(&mut buf, &pruned_compression_header, version)?;
    let compression_header_block = build_compression_header_block(&buf)?;

    let compression_header_size = compression_header_block.size(version)?;
    let mut container_size = compression_header_size;
    let mut blocks = vec![compression_header_block];
    let mut landmarks = Vec::with_capacity(slices.len());
    let mut slice_offset = compression_header_size;

    for slice in slices {
        buf.clear();

        slice::write_header(&mut buf, &slice.header, version)?;
        let slice_header_block = build_slice_header_block(&buf)?;

        landmarks.push(slice_offset);

        let mut slice_size = slice_header_block.size(version)?;
        blocks.push(slice_header_block);

        slice_size += slice.core_data_block.size(version)?;
        blocks.push(slice.core_data_block);

        for block in &slice.external_data_blocks {
            slice_size += block.size(version)?;
        }

        blocks.extend(slice.external_data_blocks);

        slice_offset += slice_size;
        container_size += slice_size;
    }

    let header = Header {
        reference_sequence_context,
        record_count,
        record_counter,
        base_count,
        block_count: blocks.len(),
        landmarks,
    };

    Ok((header, container_size, blocks))
}

fn get_container_reference_sequence_context(
    slices: &[Slice],
) -> io::Result<ReferenceSequenceContext> {
    assert!(!slices.is_empty());

    let first_slice = slices.first().expect("slices cannot be empty");
    let mut container_reference_sequence_context = first_slice.header.reference_sequence_context();

    for slice in slices.iter().skip(1) {
        let slice_reference_sequence_context = slice.header.reference_sequence_context();

        match (
            container_reference_sequence_context,
            slice_reference_sequence_context,
        ) {
            (
                ReferenceSequenceContext::Some(container_context),
                ReferenceSequenceContext::Some(slice_context),
            ) if container_context.reference_sequence_id()
                == slice_context.reference_sequence_id() =>
            {
                let alignment_start = cmp::min(
                    container_context.alignment_start(),
                    slice_context.alignment_start(),
                );

                let alignment_end = cmp::max(
                    container_context.alignment_end(),
                    slice_context.alignment_end(),
                );

                container_reference_sequence_context = ReferenceSequenceContext::some(
                    container_context.reference_sequence_id(),
                    alignment_start,
                    alignment_end,
                );
            }
            (ReferenceSequenceContext::None, ReferenceSequenceContext::None) => {}
            (ReferenceSequenceContext::Many, ReferenceSequenceContext::Many) => {}
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "invalid slice reference sequence context: expected {container_reference_sequence_context:?}, got {slice_reference_sequence_context:?}"
                    ),
                ));
            }
        }
    }

    Ok(container_reference_sequence_context)
}

fn build_compression_header_block(src: &[u8]) -> io::Result<Block> {
    Block::encode(ContentType::CompressionHeader, 0, None, src)
}

fn calculate_base_count(records: &[Record]) -> io::Result<u64> {
    let n: usize = records.iter().map(|record| record.read_length).sum();
    u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

fn build_slice_header_block(src: &[u8]) -> io::Result<Block> {
    Block::encode(ContentType::SliceHeader, 0, None, src)
}

// § 9 "End of file container" (2022-04-12)
pub static EOF_V3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

pub fn build_eof_container(version: Version) -> io::Result<Vec<u8>> {
    if version == Version::V3_0 || version == Version::V3_1 {
        return Ok(EOF_V3.to_vec());
    }

    // Build EOF dynamically for CRAM 2.x (no CRC32) and 4.0 (VLQ encoding).
    // NOTE: ref_seq_id in the EOF uses zigzag/signed encoding (matching htslib's
    // hardcoded EOF bytes), while regular container headers use unsigned VLQ/ITF8.
    use crate::io::writer::num::{
        write_i32_le, write_int, write_long, write_position, write_signed_int, write_u32_le,
    };
    use flate2::CrcWriter;

    let has_crc32 = version.has_crc32();

    // Build the empty compression header block body (without CRC32)
    let mut block_body = Vec::new();
    // compression method = none (0)
    block_body.write_all(&[0x00])?;
    // content type = compression header (1)
    block_body.write_all(&[0x01])?;
    // content id = 0
    write_int(&mut block_body, version, 0)?;
    // compressed size = 6
    write_int(&mut block_body, version, 6)?;
    // uncompressed size = 6
    write_int(&mut block_body, version, 6)?;
    // Empty compression header body: 3 components, each is a length-prefixed
    // map with 0 entries.
    block_body.write_all(&[0x01, 0x00, 0x01, 0x00, 0x01, 0x00])?;

    // Append CRC32 to block if required
    let mut block_buf = Vec::new();
    if has_crc32 {
        let mut crc_writer = CrcWriter::new(&mut block_buf);
        crc_writer.write_all(&block_body)?;
        let crc32 = crc_writer.crc().sum();
        write_u32_le(crc_writer.get_mut(), crc32)?;
    } else {
        block_buf = block_body;
    }

    // Build the container header body (without CRC32)
    let mut header_body = Vec::new();
    let block_len = i32::try_from(block_buf.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    if version.uses_vlq() {
        write_int(&mut header_body, version, block_len)?;
    } else {
        write_i32_le(&mut header_body, block_len)?;
    }
    // reference sequence id = -1 (unmapped): signed encoding in EOF
    write_signed_int(&mut header_body, version, -1)?;
    // alignment start = 0x454f46 ("EOF" magic)
    write_position(&mut header_body, version, 0x454f46)?;
    // alignment span = 0
    write_position(&mut header_body, version, 0)?;
    // record count = 0
    write_int(&mut header_body, version, 0)?;
    // record counter = 0
    write_long(&mut header_body, version, 0)?;
    // base count = 0
    write_long(&mut header_body, version, 0)?;
    // block count = 1
    write_int(&mut header_body, version, 1)?;
    // landmarks = []
    write_int(&mut header_body, version, 0)?;

    // Append CRC32 to header if required
    let mut result = Vec::new();
    if has_crc32 {
        let mut crc_writer = CrcWriter::new(&mut result);
        crc_writer.write_all(&header_body)?;
        let crc32 = crc_writer.crc().sum();
        write_u32_le(crc_writer.get_mut(), crc32)?;
    } else {
        result = header_body;
    }

    result.extend(block_buf);
    Ok(result)
}

pub fn write_eof_container<W>(writer: &mut W, version: Version) -> io::Result<()>
where
    W: Write,
{
    let eof = build_eof_container(version)?;
    writer.write_all(&eof)
}
