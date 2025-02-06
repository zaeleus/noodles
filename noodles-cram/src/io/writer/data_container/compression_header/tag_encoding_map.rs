use std::io::{self, Write};

use super::write_encoding_for_byte_array_codec;
use crate::{data_container::compression_header::TagEncodingMap, io::writer::num::write_itf8};

pub fn write_tag_encoding_map<W>(
    writer: &mut W,
    tag_encoding_map: &TagEncodingMap,
) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    let map_len = i32::try_from(tag_encoding_map.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut buf, map_len)?;

    for (&block_content_id, encoding) in tag_encoding_map.iter() {
        write_itf8(&mut buf, block_content_id)?;
        write_encoding_for_byte_array_codec(&mut buf, encoding)?;
    }

    let data_len =
        i32::try_from(buf.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, data_len)?;

    writer.write_all(&buf)?;

    Ok(())
}
