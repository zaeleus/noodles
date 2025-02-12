use std::{
    collections::HashSet,
    io::{self, Write},
};

use super::write_encoding_for_byte_array_codec;
use crate::{
    container::{
        block,
        compression_header::{
            encoding::codec::{Byte, ByteArray, Integer},
            preservation_map::tag_sets::Key,
            Encoding, TagEncodings,
        },
    },
    io::writer::{num::write_itf8, Record},
};

pub fn write_tag_encodings<W>(writer: &mut W, tag_encodings: &TagEncodings) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    let map_len = i32::try_from(tag_encodings.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut buf, map_len)?;

    for (&block_content_id, encoding) in tag_encodings.iter() {
        write_itf8(&mut buf, block_content_id)?;
        write_encoding_for_byte_array_codec(&mut buf, encoding)?;
    }

    let data_len =
        i32::try_from(buf.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, data_len)?;

    writer.write_all(&buf)?;

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
                Encoding::new(ByteArray::ByteArrayLen {
                    len_encoding: Encoding::new(Integer::External { block_content_id }),
                    value_encoding: Encoding::new(Byte::External { block_content_id }),
                }),
            )
        })
        .collect()
}
