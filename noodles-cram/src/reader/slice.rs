use std::io::{self, Read};

use crate::{
    num::{read_itf8, Itf8},
    slice,
};

pub fn read_header<R>(reader: &mut R) -> io::Result<slice::Header>
where
    R: Read,
{
    let reference_sequence_id = read_itf8(reader)?;
    let alignment_start = read_itf8(reader)?;
    let alignment_span = read_itf8(reader)?;
    let n_records = read_itf8(reader)?;
    let record_counter = read_itf8(reader)?;
    let n_blocks = read_itf8(reader)?;
    let block_content_ids = read_block_content_ids(reader)?;
    let embedded_reference_bases_block_content_id = read_itf8(reader)?;
    let reference_md5 = read_reference_md5(reader)?;
    let optional_tags = read_optional_tags(reader)?;

    Ok(slice::Header::new(
        reference_sequence_id,
        alignment_start,
        alignment_span,
        n_records,
        record_counter,
        n_blocks,
        block_content_ids,
        embedded_reference_bases_block_content_id,
        reference_md5,
        optional_tags,
    ))
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
    let len = read_itf8(reader).map(|i| i as usize)?;
    let mut buf = vec![0; len];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}
