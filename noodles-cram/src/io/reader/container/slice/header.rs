use std::{io, num::NonZeroUsize};

use bytes::{Buf, Bytes};
use noodles_core::Position;

use crate::{
    container::{
        block::{self, ContentType},
        slice::Header,
        ReferenceSequenceContext,
    },
    io::reader::{
        container::read_block,
        num::{get_itf8, get_ltf8},
    },
};

pub(super) fn get_header(src: &mut Bytes) -> io::Result<Header> {
    let block = read_block(src)?;

    if block.content_type != ContentType::SliceHeader {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block content type: expected {:?}, got {:?}",
                ContentType::SliceHeader,
                block.content_type
            ),
        ));
    }

    let mut buf = block.decode()?;
    get_header_inner(&mut buf)
}

fn get_header_inner<B>(src: &mut B) -> io::Result<Header>
where
    B: Buf,
{
    let reference_sequence_context = get_reference_sequence_context(src)?;

    let record_count = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let record_counter = get_ltf8(src).and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_count = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_content_ids = get_block_content_ids(src)?;

    let embedded_reference_bases_block_content_id =
        get_embedded_reference_bases_block_content_id(src)?;

    let reference_md5 = get_reference_md5(src)?;
    let optional_tags = get_optional_tags(src);

    let mut builder = Header::builder()
        .set_reference_sequence_context(reference_sequence_context)
        .set_record_count(record_count)
        .set_record_counter(record_counter)
        .set_block_count(block_count)
        .set_block_content_ids(block_content_ids)
        .set_reference_md5(reference_md5)
        .set_optional_tags(optional_tags);

    if let Some(id) = embedded_reference_bases_block_content_id {
        builder = builder.set_embedded_reference_bases_block_content_id(id);
    }

    Ok(builder.build())
}

fn get_reference_sequence_context<B>(src: &mut B) -> io::Result<ReferenceSequenceContext>
where
    B: Buf,
{
    const UNMAPPED: i32 = -1;
    const MULTIREF: i32 = -2;

    match get_itf8(src)? {
        UNMAPPED => {
            // Discard alignment start and span.
            get_itf8(src)?;
            get_itf8(src)?;
            Ok(ReferenceSequenceContext::None)
        }
        MULTIREF => {
            // Discard alignment start and span.
            get_itf8(src)?;
            get_itf8(src)?;
            Ok(ReferenceSequenceContext::Many)
        }
        n => {
            let reference_sequence_id =
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let alignment_start = get_itf8(src).and_then(|m| {
                usize::try_from(m)
                    .and_then(Position::try_from)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

            let alignment_span = get_itf8(src).and_then(|m| {
                usize::try_from(m)
                    .and_then(NonZeroUsize::try_from)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

            let alignment_end = alignment_start
                .checked_add(usize::from(alignment_span) - 1)
                .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;

            Ok(ReferenceSequenceContext::some(
                reference_sequence_id,
                alignment_start,
                alignment_end,
            ))
        }
    }
}

fn get_block_content_ids<B>(src: &mut B) -> io::Result<Vec<block::ContentId>>
where
    B: Buf,
{
    let len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = Vec::with_capacity(len);

    for _ in 0..len {
        let value = get_itf8(src)?;
        buf.push(value);
    }

    Ok(buf)
}

fn get_embedded_reference_bases_block_content_id<B>(
    src: &mut B,
) -> io::Result<Option<block::ContentId>>
where
    B: Buf,
{
    get_itf8(src).map(|n| match n {
        -1 => None,
        _ => Some(n),
    })
}

fn get_reference_md5<B>(src: &mut B) -> io::Result<[u8; 16]>
where
    B: Buf,
{
    let mut buf = [0; 16];

    if src.remaining() < buf.len() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    src.copy_to_slice(&mut buf);

    Ok(buf)
}

fn get_optional_tags<B>(src: &mut B) -> Vec<u8>
where
    B: Buf,
{
    let mut buf = vec![0; src.remaining()];
    src.copy_to_slice(&mut buf);
    buf
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_header_inner() -> Result<(), Box<dyn std::error::Error>> {
        let src = [
            0x02, // reference sequence ID = 2
            0x03, // alignment start = 3
            0x05, // alignment span = 5
            0x08, // number of records = 8
            0x0d, // record counter = 13
            0x01, // number of blocks = 1
            0x01, // block content ID count = 1
            0x15, // block content IDs[0] = 21
            0xff, 0xff, 0xff, 0xff, 0x0f, // embedded reference bases block content ID = -1
            0x57, 0xb2, 0x96, 0xa3, 0x16, 0x0a, 0x2c, 0xac, 0x9c, 0x83, 0x33, 0x12, 0x6f, 0xf2,
            0x7e, 0xf7, // reference MD5 (b"ACGTA")
        ];

        let actual = get_header_inner(&mut &src[..])?;

        let expected = Header::builder()
            .set_reference_sequence_context(ReferenceSequenceContext::some(
                2,
                Position::try_from(3)?,
                Position::try_from(7)?,
            ))
            .set_record_count(8)
            .set_record_counter(13)
            .set_block_count(1)
            .set_block_content_ids(vec![21])
            .set_reference_md5([
                0x57, 0xb2, 0x96, 0xa3, 0x16, 0x0a, 0x2c, 0xac, 0x9c, 0x83, 0x33, 0x12, 0x6f, 0xf2,
                0x7e, 0xf7,
            ])
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
