use std::io::{self, Write};

use crate::{
    container::compression_header::TagEncodingMap,
    num::{write_itf8, Itf8},
    writer::encoding::write_encoding,
};

#[allow(dead_code)]
pub fn write_tag_encoding_map<W>(
    writer: &mut W,
    tag_encoding_map: &TagEncodingMap,
) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    // FIXME: usize => Itf8 cast
    let map_len = tag_encoding_map.len() as Itf8;
    write_itf8(&mut buf, map_len)?;

    for (&key, encoding) in tag_encoding_map.iter() {
        write_itf8(&mut buf, key)?;
        write_encoding(&mut buf, encoding)?;
    }

    // FIXME: usize => Itf8 cast
    let data_len = buf.len() as Itf8;
    write_itf8(writer, data_len)?;

    writer.write_all(&buf)?;

    Ok(())
}
