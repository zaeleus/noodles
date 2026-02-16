use std::io;

use crate::{
    container::{
        ReferenceSequenceContext,
        block::{self, ContentType},
        slice::Header,
    },
    file_definition::Version,
    io::reader::{
        container::read_block_as,
        num::{
            read_header_int, read_long_as, read_position, read_signed_int, read_unsigned_int,
            read_unsigned_int_as,
        },
    },
};

pub(super) fn read_header(src: &mut &[u8], version: Version) -> io::Result<Header> {
    let block = read_block_as(src, ContentType::SliceHeader, version)?;
    let buf = block.decode()?;
    read_header_inner(&mut &buf[..], version)
}

fn read_header_inner(src: &mut &[u8], version: Version) -> io::Result<Header> {
    // Slice header ref_seq_id uses signed encoding (zigzag for CRAM 4.0),
    // unlike container header which uses unsigned VLQ with bitcast.
    // This matches htslib's cram_write_slice_hdr which uses varint_put32s.
    let reference_sequence_id = read_signed_int(src, version)?;

    let alignment_start = read_position(src, version)?;
    let alignment_span = read_position(src, version)?;

    let reference_sequence_context = ReferenceSequenceContext::try_from((
        reference_sequence_id,
        alignment_start,
        alignment_span,
    ))?;

    let record_count = read_unsigned_int_as(src, version)?;

    // CRAM 2.x uses ITF8 (32-bit) for record_counter per htslib;
    // CRAM 3.x+ uses LTF8 (64-bit); CRAM 4.0 uses uint7_64.
    let record_counter = if version.major() >= 3 {
        read_long_as(src, version)?
    } else {
        read_unsigned_int_as(src, version)?
    };
    let block_count = read_unsigned_int_as(src, version)?;

    let block_content_ids = read_block_content_ids(src, version)?;
    let embedded_reference_bases_block_content_id =
        read_embedded_reference_bases_block_content_id(src, version)?;
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

fn read_block_content_ids(src: &mut &[u8], version: Version) -> io::Result<Vec<block::ContentId>> {
    let len: usize = read_unsigned_int_as(src, version)?;
    (0..len).map(|_| read_unsigned_int(src, version)).collect()
}

fn read_embedded_reference_bases_block_content_id(
    src: &mut &[u8],
    version: Version,
) -> io::Result<Option<block::ContentId>> {
    // § 8.5 "Slice header block" (2024-09-04): "-1 for none".
    const MISSING: i32 = -1;

    read_header_int(src, version).map(|n| match n {
        MISSING => None,
        _ => Some(n),
    })
}

fn read_reference_md5(src: &mut &[u8]) -> io::Result<Option<[u8; 16]>> {
    let Some((buf, rest)) = src.split_first_chunk() else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    if buf.iter().all(|&b| b == 0) {
        Ok(None)
    } else {
        Ok(Some(*buf))
    }
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

        let actual = read_header_inner(&mut &src[..], Version::V3_0)?;

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
            reference_md5: Some([
                0x57, 0xb2, 0x96, 0xa3, 0x16, 0x0a, 0x2c, 0xac, 0x9c, 0x83, 0x33, 0x12, 0x6f, 0xf2,
                0x7e, 0xf7,
            ]),
            optional_tags: Vec::new(),
        };

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_inner_v2() -> Result<(), Box<dyn std::error::Error>> {
        // CRAM 2.x uses ITF8 for record_counter (32-bit) instead of LTF8.
        // The raw bytes here are identical to the V3 test because record_counter=13
        // fits in a single byte for both ITF8 and LTF8 encodings.
        let src = [
            0x02, // reference sequence ID = 2
            0x03, // alignment start = 3
            0x05, // alignment span = 5
            0x08, // number of records = 8
            0x0d, // record counter = 13 (ITF8, not LTF8)
            0x01, // number of blocks = 1
            0x01, // block content ID count = 1
            0x15, // block content IDs[0] = 21
            0xff, 0xff, 0xff, 0xff, 0x0f, // embedded reference bases block content ID = -1
            0x57, 0xb2, 0x96, 0xa3, 0x16, 0x0a, 0x2c, 0xac, 0x9c, 0x83, 0x33, 0x12, 0x6f, 0xf2,
            0x7e, 0xf7, // reference MD5 (b"ACGTA")
        ];

        let actual = read_header_inner(&mut &src[..], Version::V2_1)?;

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
            reference_md5: Some([
                0x57, 0xb2, 0x96, 0xa3, 0x16, 0x0a, 0x2c, 0xac, 0x9c, 0x83, 0x33, 0x12, 0x6f, 0xf2,
                0x7e, 0xf7,
            ]),
            optional_tags: Vec::new(),
        };

        assert_eq!(actual, expected);

        Ok(())
    }
}
