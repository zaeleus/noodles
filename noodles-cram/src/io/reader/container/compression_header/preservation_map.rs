mod substitution_matrix;
mod tag_sets;

use std::io;

use bytes::{Buf, Bytes};

use self::{substitution_matrix::get_substitution_matrix, tag_sets::get_tag_sets};
use crate::{
    container::compression_header::{preservation_map::Key, PreservationMap},
    io::reader::num::get_itf8,
};

pub(super) fn get_preservation_map(src: &mut Bytes) -> io::Result<PreservationMap> {
    let data_len = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < data_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut buf = src.split_to(data_len);
    let len = get_itf8(&mut buf)?;

    get_preservation_map_inner(&mut buf, len)
}

fn get_preservation_map_inner(src: &mut Bytes, len: i32) -> io::Result<PreservationMap> {
    let mut read_names_included = true;
    let mut ap_data_series_delta = true;
    let mut reference_required = true;
    let mut substitution_matrix = None;
    let mut tag_sets = None;

    for _ in 0..len {
        let key = get_key(src)?;

        match key {
            Key::ReadNamesIncluded => read_names_included = get_bool(src)?,
            Key::ApDataSeriesDelta => ap_data_series_delta = get_bool(src)?,
            Key::ReferenceRequired => reference_required = get_bool(src)?,
            Key::SubstitutionMatrix => {
                substitution_matrix = get_substitution_matrix(src).map(Some)?
            }
            Key::TagSets => tag_sets = get_tag_sets(src).map(Some)?,
        }
    }

    Ok(PreservationMap::new(
        read_names_included,
        ap_data_series_delta,
        reference_required,
        substitution_matrix.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing substitution matrix")
        })?,
        tag_sets.ok_or_else(|| {
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
    let n = src
        .try_get_u8()
        .map_err(|e| io::Error::new(io::ErrorKind::UnexpectedEof, e))?;

    match n {
        0 => Ok(false),
        1 => Ok(true),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid bool value",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_preservation_map() -> io::Result<()> {
        use noodles_sam::alignment::record::data::field::{Tag, Type};

        use crate::container::compression_header::preservation_map::{
            tag_sets, SubstitutionMatrix,
        };

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
            vec![vec![tag_sets::Key::new(Tag::COMMENT, Type::String)]],
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
}
