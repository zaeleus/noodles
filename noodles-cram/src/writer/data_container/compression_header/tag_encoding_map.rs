use std::io::{self, Write};

use super::write_encoding;
use crate::{
    data_container::compression_header::TagEncodingMap, num::Itf8, writer::num::write_itf8,
};

pub fn write_tag_encoding_map<W>(
    writer: &mut W,
    tag_encoding_map: &TagEncodingMap,
) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    let map_len = tag_encoding_map.len() as Itf8;
    write_itf8(&mut buf, map_len)?;

    for (&key, encoding) in tag_encoding_map.iter() {
        write_itf8(&mut buf, key)?;
        write_encoding(&mut buf, encoding)?;
    }

    let data_len = buf.len() as Itf8;
    write_itf8(writer, data_len)?;

    writer.write_all(&buf)?;

    Ok(())
}
