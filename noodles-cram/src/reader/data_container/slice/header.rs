use std::io::{self, Read};

use noodles_sam as sam;

use crate::{
    container::ReferenceSequenceId,
    data_container::slice::{self, header::EmbeddedReferenceBasesBlockContentId},
    num::Itf8,
    reader::num::{read_itf8, read_ltf8},
};

pub fn read_header<R>(reader: &mut R) -> io::Result<slice::Header>
where
    R: Read,
{
    let reference_sequence_id = read_itf8(reader).and_then(|n| {
        ReferenceSequenceId::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let alignment_start = read_itf8(reader).and_then(|n| {
        if n == 0 {
            Ok(None)
        } else {
            sam::record::Position::try_from(n)
                .map(Some)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        }
    })?;

    let alignment_span = read_itf8(reader)?;

    let record_count = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let record_counter = read_ltf8(reader)?;

    let block_count = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_content_ids = read_block_content_ids(reader)?;
    let embedded_reference_bases_block_content_id =
        read_itf8(reader).map(EmbeddedReferenceBasesBlockContentId::from)?;
    let reference_md5 = read_reference_md5(reader)?;
    let optional_tags = read_optional_tags(reader)?;

    let mut builder = slice::Header::builder()
        .set_reference_sequence_id(reference_sequence_id)
        .set_alignment_span(alignment_span)
        .set_record_count(record_count)
        .set_record_counter(record_counter)
        .set_block_count(block_count)
        .set_block_content_ids(block_content_ids)
        .set_embedded_reference_bases_block_content_id(embedded_reference_bases_block_content_id)
        .set_reference_md5(reference_md5)
        .set_optional_tags(optional_tags);

    if let Some(position) = alignment_start {
        builder = builder.set_alignment_start(position);
    }

    Ok(builder.build())
}

fn read_block_content_ids<R>(reader: &mut R) -> io::Result<Vec<Itf8>>
where
    R: Read,
{
    let len = read_itf8(reader).map(|i| i as usize)?;
    let mut buf = Vec::with_capacity(len);

    for _ in 0..len {
        let value = read_itf8(reader)?;
        buf.push(value);
    }

    Ok(buf)
}

fn read_reference_md5<R>(reader: &mut R) -> io::Result<[u8; 16]>
where
    R: Read,
{
    let mut buf = [0; 16];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

fn read_optional_tags<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_header() -> Result<(), Box<dyn std::error::Error>> {
        let data = [
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
        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = slice::Header::builder()
            .set_reference_sequence_id(ReferenceSequenceId::try_from(2)?)
            .set_alignment_start(sam::record::Position::try_from(3)?)
            .set_alignment_span(5)
            .set_record_count(8)
            .set_record_counter(13)
            .set_block_count(1)
            .set_block_content_ids(vec![21])
            .set_embedded_reference_bases_block_content_id(
                EmbeddedReferenceBasesBlockContentId::try_from(-1)?,
            )
            .set_reference_md5([
                0x57, 0xb2, 0x96, 0xa3, 0x16, 0x0a, 0x2c, 0xac, 0x9c, 0x83, 0x33, 0x12, 0x6f, 0xf2,
                0x7e, 0xf7,
            ])
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
