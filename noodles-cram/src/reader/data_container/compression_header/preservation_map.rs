use std::io::{self, BufRead, Read};

use byteorder::ReadBytesExt;
use noodles_sam::record::data::field::{value::Type, Tag};

use crate::{
    data_container::compression_header::{
        preservation_map::Key, PreservationMap, SubstitutionMatrix, TagIdsDictionary,
    },
    reader::num::read_itf8,
};

pub fn read_preservation_map<R>(reader: &mut R) -> io::Result<PreservationMap>
where
    R: Read,
{
    let data_len = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut buf = vec![0; data_len];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];
    let map_len = read_itf8(&mut buf_reader)?;

    let mut read_names_included = true;
    let mut ap_data_series_delta = true;
    let mut reference_required = true;
    let mut substitution_matrix = None;
    let mut tag_ids_dictionary = None;

    let mut key_buf = [0; 2];

    for _ in 0..map_len {
        buf_reader.read_exact(&mut key_buf)?;

        let key =
            Key::try_from(key_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        match key {
            Key::ReadNamesIncluded => {
                read_names_included = read_bool(&mut buf_reader)?;
            }
            Key::ApDataSeriesDelta => {
                ap_data_series_delta = read_bool(&mut buf_reader)?;
            }
            Key::ReferenceRequired => {
                reference_required = read_bool(&mut buf_reader)?;
            }
            Key::SubstitutionMatrix => {
                substitution_matrix = read_substitution_matrix(&mut buf_reader).map(Some)?;
            }
            Key::TagIdsDictionary => {
                tag_ids_dictionary = read_tag_ids_dictionary(&mut buf_reader).map(Some)?;
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

fn read_bool<R>(reader: &mut R) -> io::Result<bool>
where
    R: Read,
{
    match reader.read_u8() {
        Ok(0) => Ok(false),
        Ok(1) => Ok(true),
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid bool value",
        )),
        Err(e) => Err(e),
    }
}

fn read_substitution_matrix<R>(reader: &mut R) -> io::Result<SubstitutionMatrix>
where
    R: Read,
{
    let mut buf = [0; 5];
    reader.read_exact(&mut buf[..])?;
    SubstitutionMatrix::try_from(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_tag_ids_dictionary<R>(reader: &mut R) -> io::Result<TagIdsDictionary>
where
    R: Read,
{
    use crate::record::tag::Key;

    let data_len = read_itf8(reader)?;
    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];

    let mut dictionary = Vec::new();
    let mut keys_buf = Vec::new();

    loop {
        keys_buf.clear();

        match buf_reader.read_until(0x00, &mut keys_buf) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(e),
        }

        let mut line = Vec::new();

        for chunk in keys_buf.chunks_exact(3) {
            let (t0, t1, ty) = (chunk[0], chunk[1], chunk[2]);

            let tag = Tag::try_from([t0, t1])
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let ty =
                Type::try_from(ty).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let key = Key::new(tag, ty);
            line.push(key);
        }

        dictionary.push(line);
    }

    Ok(TagIdsDictionary::from(dictionary))
}

#[cfg(test)]
mod tests {
    use crate::record::tag::Key;

    use super::*;

    #[test]
    fn test_read_preservation_map() -> io::Result<()> {
        let data = [
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
        let mut reader = &data[..];
        let actual = read_preservation_map(&mut reader)?;

        let expected = PreservationMap::new(
            false,
            false,
            false,
            SubstitutionMatrix::default(),
            TagIdsDictionary::from(vec![vec![Key::new(Tag::Comment, Type::String)]]),
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_preservation_map_with_no_substitution_matrix() {
        let data = [
            0x08, // data.len = 8
            0x01, // map.len = 1
            0x54, 0x44, // key = "TD"
            0x04, 0x43, 0x4f, 0x5a, 0x00, // tag IDs dictionary = [[CO:Z]]
        ];
        let mut reader = &data[..];
        assert!(read_preservation_map(&mut reader).is_err());
    }

    #[test]
    fn test_read_preservation_map_with_no_tag_ids_dictionary() {
        let data = [
            0x08, // data.len = 8
            0x01, // map.len = 1
            0x53, 0x4d, // key = "SM"
            // [[C, G, T, N], [A, G, T, N], [A, C, T, N], [A, C, G, N], [A, C, G, T]]
            0x1b, 0x1b, 0x1b, 0x1b, 0x1b, // substitution matrix
        ];
        let mut reader = &data[..];
        assert!(read_preservation_map(&mut reader).is_err());
    }

    #[test]
    fn test_read_bool() -> io::Result<()> {
        let data = [0x00];
        let mut reader = &data[..];
        assert!(!read_bool(&mut reader)?);

        let data = [0x01];
        let mut reader = &data[..];
        assert!(read_bool(&mut reader)?);

        let data = [0x02];
        let mut reader = &data[..];
        assert!(matches!(
            read_bool(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
