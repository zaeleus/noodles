use std::{
    collections::HashSet,
    io::{self, Write},
};

use super::write_byte_array_encoding;
use crate::{
    container::{
        block,
        compression_header::{
            encoding::codec::{Byte, ByteArray, Integer},
            preservation_map::tag_sets::Key,
            Encoding, TagEncodings,
        },
    },
    io::writer::{collections::write_array, num::write_itf8, Record},
};

pub fn write_tag_encodings<W>(writer: &mut W, tag_encodings: &TagEncodings) -> io::Result<()>
where
    W: Write,
{
    let buf = encode(tag_encodings)?;
    write_array(writer, &buf)
}

fn encode(tag_encodings: &TagEncodings) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    encode_inner(&mut buf, tag_encodings)?;
    Ok(buf)
}

fn encode_inner<W>(writer: &mut W, tag_encodings: &TagEncodings) -> io::Result<()>
where
    W: Write,
{
    let len = i32::try_from(tag_encodings.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, len)?;

    for (&block_content_id, encoding) in tag_encodings.iter() {
        write_itf8(writer, block_content_id)?;
        write_byte_array_encoding(writer, encoding)?;
    }

    Ok(())
}

pub(super) fn build_tag_encodings(records: &[Record]) -> TagEncodings {
    let mut block_content_ids = HashSet::new();

    for record in records {
        for (tag, value) in &record.data {
            let key = Key::new(*tag, value.ty());
            let block_content_id = block::ContentId::from(key);
            block_content_ids.insert(block_content_id);
        }
    }

    block_content_ids
        .into_iter()
        .map(|block_content_id| {
            (
                block_content_id,
                Encoding::new(ByteArray::ByteArrayLength {
                    len_encoding: Encoding::new(Integer::External { block_content_id }),
                    value_encoding: Encoding::new(Byte::External { block_content_id }),
                }),
            )
        })
        .collect()
}
