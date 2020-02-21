use std::{
    convert::TryFrom,
    io::{self, Read},
};

use crate::{
    block::ContentType,
    num::{read_itf8, Itf8},
    slice, Block, Slice,
};

use super::block::read_block;

pub fn read_header<R>(reader: &mut R) -> io::Result<slice::Header>
where
    R: Read,
{
    let mut block = Block::default();
    read_block(reader, &mut block)?;

    let data = block.decompressed_data();
    let mut data_reader = &data[..];

    let reference_sequence_id = read_itf8(&mut data_reader)?;
    let alignment_start = read_itf8(&mut data_reader)?;
    let alignment_span = read_itf8(&mut data_reader)?;
    let n_records = read_itf8(&mut data_reader)?;
    let record_counter = read_itf8(&mut data_reader)?;
    let n_blocks = read_itf8(&mut data_reader)?;
    let block_content_ids = read_block_content_ids(&mut data_reader)?;
    let embedded_reference_bases_block_content_id = read_itf8(&mut data_reader)?;
    let reference_md5 = read_reference_md5(&mut data_reader)?;
    let optional_tags = read_optional_tags(&mut data_reader)?;

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

pub fn read_blocks<R>(reader: &mut R, slice: &mut Slice) -> io::Result<()>
where
    R: Read,
{
    for _ in 0..slice.header().n_blocks() {
        let mut block = Block::default();
        read_block(reader, &mut block)?;

        let content_type =
            ContentType::try_from(block.content_type()).expect("invalid content type");

        match content_type {
            ContentType::CoreData => {
                *slice.core_data_block_mut() = block;
            }
            ContentType::ExternalData => {
                slice.add_external_block(block);
            }
            ty => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "invalid block content type: expected CoreData or ExternalData, got {:?}",
                        ty
                    ),
                ))
            }
        }
    }

    Ok(())
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
        Err(e) => match e.kind() {
            io::ErrorKind::UnexpectedEof => return Ok(Vec::new()),
            _ => return Err(e),
        },
    };

    let mut buf = vec![0; len];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}
