use std::io;

use bytes::{Buf, Bytes};
use noodles_sam::alignment::record::data::field::{Tag, Type};

use crate::{
    container::compression_header::preservation_map::{tag_sets::Key, TagSets},
    io::reader::num::get_itf8,
};

pub(super) fn get_tag_sets(src: &mut Bytes) -> io::Result<TagSets> {
    const NUL: u8 = 0x00;

    let data_len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < data_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut buf = src.split_to(data_len);

    let mut sets = Vec::new();

    while let Some(i) = buf.iter().position(|&b| b == NUL) {
        let keys_buf = buf.split_to(i);
        buf.advance(1); // Discard the NUL terminator.

        let mut line = Vec::new();

        for chunk in keys_buf.chunks_exact(3) {
            let (t0, t1, ty) = (chunk[0], chunk[1], chunk[2]);

            let tag = Tag::new(t0, t1);
            let ty = get_type(ty)?;
            let key = Key::new(tag, ty);

            line.push(key);
        }

        sets.push(line);
    }

    Ok(sets)
}

fn get_type(n: u8) -> io::Result<Type> {
    match n {
        b'A' => Ok(Type::Character),
        b'c' => Ok(Type::Int8),
        b'C' => Ok(Type::UInt8),
        b's' => Ok(Type::Int16),
        b'S' => Ok(Type::UInt16),
        b'i' => Ok(Type::Int32),
        b'I' => Ok(Type::UInt32),
        b'f' => Ok(Type::Float),
        b'Z' => Ok(Type::String),
        b'H' => Ok(Type::Hex),
        b'B' => Ok(Type::Array),
        _ => Err(io::Error::from(io::ErrorKind::InvalidData)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_type() -> io::Result<()> {
        assert_eq!(get_type(b'A')?, Type::Character);
        assert_eq!(get_type(b'c')?, Type::Int8);
        assert_eq!(get_type(b'C')?, Type::UInt8);
        assert_eq!(get_type(b's')?, Type::Int16);
        assert_eq!(get_type(b'S')?, Type::UInt16);
        assert_eq!(get_type(b'i')?, Type::Int32);
        assert_eq!(get_type(b'I')?, Type::UInt32);
        assert_eq!(get_type(b'f')?, Type::Float);
        assert_eq!(get_type(b'Z')?, Type::String);
        assert_eq!(get_type(b'H')?, Type::Hex);
        assert_eq!(get_type(b'B')?, Type::Array);

        assert!(matches!(
            get_type(b'n'),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
