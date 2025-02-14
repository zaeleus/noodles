use std::io;

use crate::{
    container::{
        block::{self, ContentType},
        slice::Header,
        ReferenceSequenceContext,
    },
    io::reader::{
        container::read_block,
        num::{read_itf8, read_itf8_as, read_ltf8_as},
        split_at_checked,
    },
};

pub(super) fn read_header(src: &mut &[u8]) -> io::Result<Header> {
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

    let buf = block.decode()?;
    read_header_inner(&mut &buf[..])
}

fn read_header_inner(src: &mut &[u8]) -> io::Result<Header> {
    let reference_sequence_id = read_itf8(src)?;
    let alignment_start = read_itf8(src)?;
    let alignment_span = read_itf8(src)?;

    let reference_sequence_context = ReferenceSequenceContext::try_from((
        reference_sequence_id,
        alignment_start,
        alignment_span,
    ))?;

    let record_count = read_itf8_as(src)?;
    let record_counter = read_ltf8_as(src)?;
    let block_count = read_itf8_as(src)?;

    let block_content_ids = read_block_content_ids(src)?;
    let embedded_reference_bases_block_content_id =
        read_embedded_reference_bases_block_content_id(src)?;
    let reference_md5 = read_reference_md5(src)?;
    let optional_tags = read_optional_tags(src);

    Ok(Header {
        reference_sequence_context,
        record_count,
        record_counter,
        block_count,
        block_content_ids,
        embedded_reference_bases_block_content_id,
        reference_md5,
        optional_tags,
    })
}

fn read_block_content_ids(src: &mut &[u8]) -> io::Result<Vec<block::ContentId>> {
    let len: usize = read_itf8_as(src)?;
    (0..len).map(|_| read_itf8(src)).collect()
}

fn read_embedded_reference_bases_block_content_id(
    src: &mut &[u8],
) -> io::Result<Option<block::ContentId>> {
    // § 8.5 "Slice header block" (2024-09-04): "-1 for none".
    const MISSING: i32 = -1;

    read_itf8(src).map(|n| match n {
        MISSING => None,
        _ => Some(n),
    })
}

fn read_reference_md5(src: &mut &[u8]) -> io::Result<[u8; 16]> {
    const SIZE: usize = 16;

    let Some((buf, rest)) = split_at_checked(src, SIZE) else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    // SAFETY: `buf.len() == 16`.
    Ok(buf.try_into().unwrap())
}

fn read_optional_tags(src: &mut &[u8]) -> Vec<u8> {
    let (buf, rest) = src.split_at(src.len());
    *src = rest;
    buf.into()
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;

    #[test]
    fn test_read_header_inner() -> Result<(), Box<dyn std::error::Error>> {
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

        let actual = read_header_inner(&mut &src[..])?;

        let expected = Header {
            reference_sequence_context: ReferenceSequenceContext::some(
                2,
                Position::try_from(3)?,
                Position::try_from(7)?,
            ),
            record_count: 8,
            record_counter: 13,
            block_count: 1,
            block_content_ids: vec![21],
            embedded_reference_bases_block_content_id: None,
            reference_md5: [
                0x57, 0xb2, 0x96, 0xa3, 0x16, 0x0a, 0x2c, 0xac, 0x9c, 0x83, 0x33, 0x12, 0x6f, 0xf2,
                0x7e, 0xf7,
            ],
            optional_tags: Vec::new(),
        };

        assert_eq!(actual, expected);

        Ok(())
    }
}
