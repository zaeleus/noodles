mod substitution_matrix;
mod tag_sets;

use std::io;

use self::{substitution_matrix::read_substitution_matrix, tag_sets::read_tag_sets};
use crate::{
    container::compression_header::{PreservationMap, preservation_map::Key},
    file_definition::Version,
    io::reader::collections::read_map,
};

pub(super) fn read_preservation_map(
    src: &mut &[u8],
    version: Version,
) -> io::Result<PreservationMap> {
    let (mut buf, len) = read_map(src, version)?;
    read_preservation_map_inner(&mut buf, len, version)
}

fn read_preservation_map_inner(
    src: &mut &[u8],
    len: usize,
    version: Version,
) -> io::Result<PreservationMap> {
    let mut records_have_names = true;
    let mut alignment_starts_are_deltas = true;
    let mut external_reference_sequence_is_required = true;
    let mut substitution_matrix = None;
    let mut tag_sets = None;
    let mut qs_seq_orient = true;

    for _ in 0..len {
        let key = read_key(src)?;

        match key {
            Key::RecordsHaveNames => records_have_names = read_bool(src)?,
            Key::AlignmentStartsAreDeltas => alignment_starts_are_deltas = read_bool(src)?,
            Key::ExternalReferenceSequenceIsRequired => {
                external_reference_sequence_is_required = read_bool(src)?
            }
            Key::SubstitutionMatrix => {
                substitution_matrix = read_substitution_matrix(src).map(Some)?
            }
            Key::TagSets => tag_sets = read_tag_sets(src, version).map(Some)?,
            Key::QualityScoreOrientation => qs_seq_orient = read_quality_score_orientation(src)?,
        }
    }

    let mut map = PreservationMap::new(
        records_have_names,
        alignment_starts_are_deltas,
        external_reference_sequence_is_required,
        substitution_matrix.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing substitution matrix")
        })?,
        tag_sets.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing tag IDs dictionary")
        })?,
    );

    map.qs_seq_orient = qs_seq_orient;

    Ok(map)
}

fn read_key(src: &mut &[u8]) -> io::Result<Key> {
    let (buf, rest) = src
        .split_first_chunk()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Key::try_from(*buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

// CRAM 4.0: Quality score orientation (QO). A single byte:
// 0 = original/sequencing orientation, 1 = alignment orientation.
// Intentionally lenient (any non-zero = true) to match htslib behavior,
// rather than strict 0x00/0x01 validation used by read_bool.
fn read_quality_score_orientation(src: &mut &[u8]) -> io::Result<bool> {
    let (&n, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(n != 0)
}

// § 2.3 "Writing bytes to a byte stream" (2024-09-04): "Boolean is written as 1-byte with 0x0
// being 'false' and 0x1 being 'true'."
fn read_bool(src: &mut &[u8]) -> io::Result<bool> {
    const FALSE: u8 = 0x00;
    const TRUE: u8 = 0x01;

    let (n, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

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
        ];

        let actual = read_preservation_map(&mut &src[..], Version::V3_0)?;

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
            0x54, 0x44, // key = "TD"
            0x04, 0x43, 0x4f, 0x5a, 0x00, // tag IDs dictionary = [[CO:Z]]
        ];

        assert!(matches!(
            read_preservation_map(&mut &src[..], Version::V3_0),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }

    #[test]
    fn test_read_preservation_map_with_no_tag_ids_dictionary() {
        let src = [
            0x08, // data.len = 8
            0x01, // map.len = 1
            0x53, 0x4d, // key = "SM"
            // [[C, G, T, N], [A, G, T, N], [A, C, T, N], [A, C, G, N], [A, C, G, T]]
            0x1b, 0x1b, 0x1b, 0x1b, 0x1b, // substitution matrix
        ];

        assert!(matches!(
            read_preservation_map(&mut &src[..], Version::V3_0),
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
