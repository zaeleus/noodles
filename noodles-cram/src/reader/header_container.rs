mod header;

use std::{
    io::{self, Read},
    str,
};

use bytes::{Buf, Bytes, BytesMut};

use self::header::read_header;

pub fn read_header_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<String>
where
    R: Read,
{
    let len = read_header(reader)?;

    buf.resize(len, 0);
    reader.read_exact(buf)?;
    let mut buf = buf.split().freeze();

    read_raw_sam_header_from_block(&mut buf)
}

pub fn read_raw_sam_header_from_block(src: &mut Bytes) -> io::Result<String> {
    use super::container::read_block;
    use crate::container::block::ContentType;

    const EXPECTED_CONTENT_TYPE: ContentType = ContentType::FileHeader;

    let block = read_block(src)?;

    if block.content_type() != EXPECTED_CONTENT_TYPE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block content type: expected {:?}, got {:?}",
                EXPECTED_CONTENT_TYPE,
                block.content_type()
            ),
        ));
    }

    let mut data = block.decompressed_data()?;
    let _header_len = data.get_i32_le();

    str::from_utf8(&data[..])
        .map(|s| s.into())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
