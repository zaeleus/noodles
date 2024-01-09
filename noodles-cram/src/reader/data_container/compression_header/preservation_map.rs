use std::io;

use bytes::{Buf, Bytes};
use noodles_sam::record::data::field::{Tag, Type};

use crate::{
    data_container::compression_header::{
        preservation_map::{tag_ids_dictionary, Key},
        PreservationMap, SubstitutionMatrix, TagIdsDictionary,
    },
    reader::num::get_itf8,
};

pub(super) fn get_preservation_map(src: &mut Bytes) -> io::Result<PreservationMap> {
    let data_len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < data_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut buf = src.split_to(data_len);

    let map_len = get_itf8(&mut buf)?;

    let mut read_names_included = true;
    let mut ap_data_series_delta = true;
    let mut reference_required = true;
    let mut substitution_matrix = None;
    let mut tag_ids_dictionary = None;

    for _ in 0..map_len {
        let key = get_key(&mut buf)?;

        match key {
            Key::ReadNamesIncluded => {
                read_names_included = get_bool(&mut buf)?;
            }
            Key::ApDataSeriesDelta => {
                ap_data_series_delta = get_bool(&mut buf)?;
            }
            Key::ReferenceRequired => {
                reference_required = get_bool(&mut buf)?;
            }
            Key::SubstitutionMatrix => {
                substitution_matrix = get_substitution_matrix(&mut buf).map(Some)?;
            }
            Key::TagIdsDictionary => {
                tag_ids_dictionary = get_tag_ids_dictionary(&mut buf).map(Some)?;
            }
        }
    }

    Ok(PreservationMap::new(
        read_names_included,
        ap_data_series_delta,
        reference_required,
        substitution_matrix.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing substitution matrix")
        })?,
        tag_ids_dictionary.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing tag IDs dictionary")
        })?,
    ))
}

fn get_key<B>(src: &mut B) -> io::Result<Key>
where
    B: Buf,
{
    let mut buf = [0; 2];

    if src.remaining() < buf.len() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    src.copy_to_slice(&mut buf);

    Key::try_from(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn get_bool<B>(src: &mut B) -> io::Result<bool>
where
    B: Buf,
{
    if !src.has_remaining() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_u8() {
        0 => Ok(false),
        1 => Ok(true),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid bool value",
        )),
    }
}

fn get_substitution_matrix<B>(src: &mut B) -> io::Result<SubstitutionMatrix>
where
    B: Buf,
{
    let mut buf = [0; 5];

    if src.remaining() < buf.len() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    src.copy_to_slice(&mut buf);

    SubstitutionMatrix::try_from(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn get_tag_ids_dictionary(src: &mut Bytes) -> io::Result<TagIdsDictionary> {
    const NUL: u8 = 0x00;

    let data_len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < data_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut buf = src.split_to(data_len);

    let mut dictionary = Vec::new();

    while let Some(i) = buf.iter().position(|&b| b == NUL) {
        let keys_buf = buf.split_to(i);
        buf.advance(1); // Discard the NUL terminator.

        let mut line = Vec::new();

        for chunk in keys_buf.chunks_exact(3) {
            let (t0, t1, ty) = (chunk[0], chunk[1], chunk[2]);

            let tag = Tag::try_from([t0, t1])
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let ty = get_type(ty)?;

            let key = tag_ids_dictionary::Key::new(tag, ty);
            line.push(key);
        }

        dictionary.push(line);
    }

    Ok(TagIdsDictionary::from(dictionary))
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
    fn test_get_preservation_map() -> io::Result<()> {
        use noodles_sam::record::data::field::tag;

        let mut data = Bytes::from_static(&[
            0x18, // data.len = 24
            0x05, // map.len = 5
            0x52, 0x4e, // key = "RN"
            0x00, // map["RN"] = false
            0x41, 0x50, // key = "AP"
            0x00, // map["AP"] = false
            0x52, 0x52, // key = "RR"
            0x00, // map["RR"] = false
            0x53, 0x4d, // key = "SM"
            // [[C, G, T, N], [A, G, T, N], [A, C, T, N], [A, C, G, N], [A, C, G, T]]
            0x1b, 0x1b, 0x1b, 0x1b, 0x1b, // substitution matrix
            0x54, 0x44, // key = "TD"
            0x04, 0x43, 0x4f, 0x5a, 0x00, // tag IDs dictionary = [[CO:Z]]
        ]);

        let actual = get_preservation_map(&mut data)?;

        let expected = PreservationMap::new(
            false,
            false,
            false,
            SubstitutionMatrix::default(),
            TagIdsDictionary::from(vec![vec![tag_ids_dictionary::Key::new(
                tag::COMMENT,
                Type::String,
            )]]),
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_get_preservation_map_with_no_substitution_matrix() {
        let mut data = Bytes::from_static(&[
            0x08, // data.len = 8
            0x01, // map.len = 1
            0x54, 0x44, // key = "TD"
            0x04, 0x43, 0x4f, 0x5a, 0x00, // tag IDs dictionary = [[CO:Z]]
        ]);

        assert!(matches!(
            get_preservation_map(&mut data),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }

    #[test]
    fn test_get_preservation_map_with_no_tag_ids_dictionary() {
        let mut data = Bytes::from_static(&[
            0x08, // data.len = 8
            0x01, // map.len = 1
            0x53, 0x4d, // key = "SM"
            // [[C, G, T, N], [A, G, T, N], [A, C, T, N], [A, C, G, N], [A, C, G, T]]
            0x1b, 0x1b, 0x1b, 0x1b, 0x1b, // substitution matrix
        ]);

        assert!(matches!(
            get_preservation_map(&mut data),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }

    #[test]
    fn test_get_bool() -> io::Result<()> {
        let data = [0x00];
        let mut reader = &data[..];
        assert!(!get_bool(&mut reader)?);

        let data = [0x01];
        let mut reader = &data[..];
        assert!(get_bool(&mut reader)?);

        let data = [0x02];
        let mut reader = &data[..];
        assert!(matches!(
            get_bool(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

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
