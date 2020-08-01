use std::{
    convert::TryFrom,
    io::{self, BufRead, Read},
};

use byteorder::ReadBytesExt;

use crate::{
    container::compression_header::{
        preservation_map::Key, PreservationMap, SubstitutionMatrix, TagIdsDictionary,
    },
    num::read_itf8,
};

pub fn read_preservation_map<R>(reader: &mut R) -> io::Result<PreservationMap>
where
    R: Read,
{
    let data_len = read_itf8(reader)?;

    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];
    let map_len = read_itf8(&mut buf_reader)?;

    let mut map = PreservationMap::default();
    let mut key_buf = [0; 2];

    for _ in 0..map_len {
        buf_reader.read_exact(&mut key_buf)?;

        let key = Key::try_from(&key_buf[..])
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        match key {
            Key::ReadNamesIncluded => {
                *map.read_names_included_mut() = read_bool(&mut buf_reader)?;
            }
            Key::ApDataSeriesDelta => {
                *map.ap_data_series_delta_mut() = read_bool(&mut buf_reader)?;
            }
            Key::ReferenceRequired => {
                *map.reference_required_mut() = read_bool(&mut buf_reader)?;
            }
            Key::SubstitutionMatrix => {
                let mut buf = [0; 5];
                buf_reader.read_exact(&mut buf[..])?;

                let matrix = SubstitutionMatrix::try_from(&buf[..])
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                *map.substitution_matrix_mut() = matrix;
            }
            Key::TagIdsDictionary => {
                *map.tag_ids_dictionary_mut() = read_tag_ids_dictionary(&mut buf_reader)?;
            }
        }
    }

    Ok(map)
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

fn read_tag_ids_dictionary<R>(reader: &mut R) -> io::Result<TagIdsDictionary>
where
    R: Read,
{
    let data_len = read_itf8(reader)?;
    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];

    let mut dictionary = TagIdsDictionary::new();
    let mut keys_buf = Vec::new();

    loop {
        keys_buf.clear();

        match buf_reader.read_until(0x00, &mut keys_buf) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(e),
        }

        let keys: Vec<_> = keys_buf.chunks_exact(3).map(|k| k.to_vec()).collect();

        dictionary.push(keys);
    }

    Ok(dictionary)
}
