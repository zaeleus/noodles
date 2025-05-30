mod block;
pub(crate) mod compression_header;
mod header;
pub(crate) mod slice;

use std::{
    cmp,
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
use super::{DEFAULT_RECORDS_PER_SLICE, Options, Record};
use crate::container::{Header, ReferenceSequenceContext, block::ContentType};

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

    let (header, container_size, blocks) = build_container(
        reference_sequence_repository,
        options,
        header,
        record_counter,
        records,
    )?;

    write_header(writer, &header, container_size)?;

    for block in blocks {
        write_block(writer, &block)?;
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
    let mut slices = Vec::new();
    let mut slice_record_counter = record_counter;

    let compression_header = build_compression_header(options, records);

    for chunk in records.chunks_mut(DEFAULT_RECORDS_PER_SLICE) {
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

    let reference_sequence_context = get_container_reference_sequence_context(&slices)?;
    let record_count = records.len();
    let base_count = calculate_base_count(records)?;

    let mut buf = Vec::new();

    write_compression_header(&mut buf, &compression_header)?;
    let compression_header_block = build_compression_header_block(&buf)?;

    let mut container_size = compression_header_block.size()?;
    let mut blocks = vec![compression_header_block];
    let mut landmarks = Vec::with_capacity(slices.len());

    for slice in slices {
        buf.clear();

        slice::write_header(&mut buf, &slice.header)?;
        let slice_header_block = build_slice_header_block(&buf)?;

        let mut slice_size = slice_header_block.size()?;
        blocks.push(slice_header_block);

        slice_size += slice.core_data_block.size()?;
        blocks.push(slice.core_data_block);

        for block in &slice.external_data_blocks {
            slice_size += block.size()?;
        }

        blocks.extend(slice.external_data_blocks);

        let last_landmark = landmarks.last().copied().unwrap_or(0);
        let landmark = last_landmark + slice_size;
        landmarks.push(landmark);

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
pub static EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

pub fn write_eof_container<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&EOF)
}
