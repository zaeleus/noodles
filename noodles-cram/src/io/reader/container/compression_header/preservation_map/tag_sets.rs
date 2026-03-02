use std::io;

use noodles_sam::alignment::record::data::field::{Tag, Type};

use crate::{
    container::compression_header::preservation_map::{TagSets, tag_sets::Key},
    file_definition::Version,
    io::reader::collections::read_array,
};

pub(super) fn read_tag_sets(src: &mut &[u8], version: Version) -> io::Result<TagSets> {
    let mut buf = read_array(src, version)?;
    read_tag_sets_inner(&mut buf)
}

fn read_tag_sets_inner(src: &mut &[u8]) -> io::Result<TagSets> {
    const NUL: u8 = 0x00;

    let mut sets = Vec::new();

    while let Some(i) = src.iter().position(|&b| b == NUL) {
        let (buf, rest) = src.split_at(i);

        *src = &rest[1..];

        let mut line = Vec::new();

        for chunk in buf.chunks_exact(3) {
            let (t0, t1, ty) = (chunk[0], chunk[1], chunk[2]);

            let tag = Tag::new(t0, t1);
            let ty = decode_type(ty)?;
            let key = Key::new(tag, ty);

            line.push(key);
        }

        sets.push(line);
    }

    Ok(sets)
}

fn decode_type(n: u8) -> io::Result<Type> {
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
    fn test_decode_type() -> io::Result<()> {
        assert_eq!(decode_type(b'A')?, Type::Character);
        assert_eq!(decode_type(b'c')?, Type::Int8);
        assert_eq!(decode_type(b'C')?, Type::UInt8);
        assert_eq!(decode_type(b's')?, Type::Int16);
        assert_eq!(decode_type(b'S')?, Type::UInt16);
        assert_eq!(decode_type(b'i')?, Type::Int32);
        assert_eq!(decode_type(b'I')?, Type::UInt32);
        assert_eq!(decode_type(b'f')?, Type::Float);
        assert_eq!(decode_type(b'Z')?, Type::String);
        assert_eq!(decode_type(b'H')?, Type::Hex);
        assert_eq!(decode_type(b'B')?, Type::Array);

        assert!(matches!(
            decode_type(b'n'),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
