mod substitution_matrix;
mod tag_sets;

use std::io;

use self::{substitution_matrix::read_substitution_matrix, tag_sets::read_tag_sets};
use crate::{container::compression_header::PreservationMap, io::reader::collections::read_map};

pub(super) fn read_preservation_map(src: &mut &[u8]) -> io::Result<PreservationMap> {
    let (mut buf, len) = read_map(src)?;
    read_preservation_map_inner(&mut buf, len)
}

fn read_preservation_map_inner(src: &mut &[u8], len: usize) -> io::Result<PreservationMap> {
    use crate::container::compression_header::preservation_map::key;

    let mut records_have_names = true;
    let mut alignment_starts_are_deltas = true;
    let mut external_reference_sequence_is_required = true;
    let mut substitution_matrix = None;
    let mut tag_sets = None;

    for _ in 0..len {
        let key = read_key(src)?;

        match key {
            key::RECORDS_HAVE_NAMES => records_have_names = read_bool(src)?,
            key::ALIGNMENT_STARTS_ARE_DELTAS => alignment_starts_are_deltas = read_bool(src)?,
            key::EXTERNAL_REFERENCE_SEQUENCE_IS_REQUIRED => {
                external_reference_sequence_is_required = read_bool(src)?
            }
            key::SUBSTITUTION_MATRIX => {
                substitution_matrix = read_substitution_matrix(src).map(Some)?
            }
            key::TAG_SETS => tag_sets = read_tag_sets(src).map(Some)?,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid preservation map key",
                ));
            }
        }
    }

    Ok(PreservationMap::new(
        records_have_names,
        alignment_starts_are_deltas,
        external_reference_sequence_is_required,
        substitution_matrix.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing substitution matrix")
        })?,
        tag_sets.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing tag IDs dictionary")
        })?,
    ))
}

fn read_key<'a>(src: &mut &'a [u8]) -> io::Result<&'a [u8; 2]> {
    let (key, rest) = src
        .split_first_chunk()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(key)
}

// § 2.3 "Writing bytes to a byte stream" (2024-09-04): "Boolean is written as 1-byte with 0x0
// being 'false' and 0x1 being 'true'."
fn read_bool(src: &mut &[u8]) -> io::Result<bool> {
    const FALSE: u8 = 0x00;
    const TRUE: u8 = 0x01;

    let n = src
        .split_off_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    match *n {
        FALSE => Ok(false),
        TRUE => Ok(true),
        _ => Err(io::Error::new(io::ErrorKind::InvalidData, "invalid bool")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_preservation_map() -> io::Result<()> {
        use noodles_sam::alignment::record::data::field::{Tag, Type};

        use crate::container::compression_header::preservation_map::{
            SubstitutionMatrix, tag_sets,
        };

        let src = [
            0x18, // data.len = 24
            0x05, // map.len = 5
            b'R', b'N', // key = "RN"
            0x00, // map["RN"] = false
            b'A', b'P', // key = "AP"
            0x00, // map["AP"] = false
            b'R', b'R', // key = "RR"
            0x00, // map["RR"] = false
            b'S', b'M', // key = "SM"
            // [[C, G, T, N], [A, G, T, N], [A, C, T, N], [A, C, G, N], [A, C, G, T]]
            0x1b, 0x1b, 0x1b, 0x1b, 0x1b, // substitution matrix
            b'T', b'D', // key = "TD"
            0x04, b'C', b'O', b'Z', 0x00, // tag IDs dictionary = [[CO:Z]]
        ];

        let actual = read_preservation_map(&mut &src[..])?;

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
    fn test_read_preservation_map_with_no_substitution_matrix() {
        let src = [
            0x08, // data.len = 8
            0x01, // map.len = 1
            b'T', b'D', // key = "TD"
            0x04, b'C', b'O', b'Z', 0x00, // tag IDs dictionary = [[CO:Z]]
        ];

        assert!(matches!(
            read_preservation_map(&mut &src[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }

    #[test]
    fn test_read_preservation_map_with_no_tag_ids_dictionary() {
        let src = [
            0x08, // data.len = 8
            0x01, // map.len = 1
            b'S', b'M', // key = "SM"
            // [[C, G, T, N], [A, G, T, N], [A, C, T, N], [A, C, G, N], [A, C, G, T]]
            0x1b, 0x1b, 0x1b, 0x1b, 0x1b, // substitution matrix
        ];

        assert!(matches!(
            read_preservation_map(&mut &src[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }

    #[test]
    fn test_read_bool() -> io::Result<()> {
        assert!(!read_bool(&mut &[0x00][..])?);
        assert!(read_bool(&mut &[0x01][..])?);

        assert!(matches!(
            read_bool(&mut &[0x02][..]),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
