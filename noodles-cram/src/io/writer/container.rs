mod block;
pub(crate) mod compression_header;
mod header;
mod slice;

use std::{
    cmp,
    io::{self, Write},
};

use self::compression_header::write_compression_header;
pub use self::{block::write_block, header::write_header};
use crate::{
    container::{Block, Header, ReferenceSequenceContext, Slice},
    Container,
};

pub fn write_container<W>(writer: &mut W, container: &Container, base_count: u64) -> io::Result<()>
where
    W: Write,
{
    let (header, container_size, blocks) = build_container(container, base_count)?;

    write_header(writer, &header, container_size)?;

    for block in blocks {
        write_block(writer, &block)?;
    }

    Ok(())
}

fn build_container(
    container: &Container,
    base_count: u64,
) -> io::Result<(Header, usize, Vec<Block>)> {
    use crate::container::block::ContentType;

    let mut buf = Vec::new();
    write_compression_header(&mut buf, container.compression_header())?;

    let block = Block::builder()
        .set_content_type(ContentType::CompressionHeader)
        .set_uncompressed_len(buf.len())
        .set_data(buf.into())
        .build();

    let mut blocks = vec![block];
    let mut landmarks = Vec::new();

    let container_reference_sequence_context =
        build_container_reference_sequence_context(container.slices())?;

    let mut container_record_count = 0;
    let container_record_counter = container
        .slices()
        .first()
        .map(|s| s.header().record_counter())
        .expect("no slices in builder");

    for slice in container.slices() {
        let slice_header = slice.header();

        container_record_count += slice_header.record_count() as i32;

        let mut slice_len = 0;

        let mut slice_header_buf = Vec::new();
        self::slice::write_header(&mut slice_header_buf, slice.header())?;

        let slice_header_block = Block::builder()
            .set_content_type(ContentType::SliceHeader)
            .set_uncompressed_len(slice_header_buf.len())
            .set_data(slice_header_buf.into())
            .build();

        slice_len += slice_header_block.len();
        blocks.push(slice_header_block);

        blocks.push(slice.core_data_block().clone());
        slice_len += slice.core_data_block().len();

        for external_block in slice.external_blocks() {
            blocks.push(external_block.clone());
            slice_len += external_block.len();
        }

        let last_landmark = landmarks.last().copied().unwrap_or(0);
        let landmark = last_landmark + slice_len;
        landmarks.push(landmark);
    }

    let len = blocks.iter().map(|b| b.len()).sum();

    let header = Header::builder()
        .set_reference_sequence_context(container_reference_sequence_context)
        .set_record_count(container_record_count)
        .set_record_counter(container_record_counter)
        .set_base_count(base_count)
        .set_block_count(blocks.len())
        .set_landmarks(landmarks)
        .build();

    Ok((header, len, blocks))
}

fn build_container_reference_sequence_context(
    slices: &[Slice],
) -> io::Result<ReferenceSequenceContext> {
    assert!(!slices.is_empty());

    let first_slice = slices.first().expect("slices cannot be empty");
    let mut container_reference_sequence_context =
        first_slice.header().reference_sequence_context();

    for slice in slices.iter().skip(1) {
        let slice_reference_sequence_context = slice.header().reference_sequence_context();

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
