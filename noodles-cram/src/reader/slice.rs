use std::{
    convert::TryFrom,
    io::{self, Read},
};

use crate::{
    container::{
        slice::{self, header::EmbeddedReferenceBasesBlockContentId},
        ReferenceSequenceId,
    },
    num::{read_itf8, read_ltf8, Itf8},
};

pub fn read_header<R>(reader: &mut R) -> io::Result<slice::Header>
where
    R: Read,
{
    let reference_sequence_id = read_itf8(reader).and_then(|n| {
        ReferenceSequenceId::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let alignment_start = read_itf8(reader)?;
    let alignment_span = read_itf8(reader)?;
    let record_count = read_itf8(reader)?;
    let record_counter = read_ltf8(reader)?;
    let block_count = read_itf8(reader)?;
    let block_content_ids = read_block_content_ids(reader)?;
    let embedded_reference_bases_block_content_id =
        read_itf8(reader).map(EmbeddedReferenceBasesBlockContentId::from)?;
    let reference_md5 = read_reference_md5(reader)?;
    let optional_tags = read_optional_tags(reader)?;

    Ok(slice::Header::builder()
        .set_reference_sequence_id(reference_sequence_id)
        .set_alignment_start(alignment_start)
        .set_alignment_span(alignment_span)
        .set_record_count(record_count)
        .set_record_counter(record_counter)
        .set_block_count(block_count)
        .set_block_content_ids(block_content_ids)
        .set_embedded_reference_bases_block_content_id(embedded_reference_bases_block_content_id)
        .set_reference_md5(reference_md5)
        .set_optional_tags(optional_tags)
        .build())
}

fn read_block_content_ids<R>(reader: &mut R) -> io::Result<Vec<Itf8>>
where
    R: Read,
{
    let len = read_itf8(reader).map(|i| i as usize)?;
    let mut buf = vec![0; len];

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
    let len = match read_itf8(reader) {
        Ok(len) => len as usize,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(Vec::new()),
        Err(e) => return Err(e),
    };

    let mut buf = vec![0; len];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}
