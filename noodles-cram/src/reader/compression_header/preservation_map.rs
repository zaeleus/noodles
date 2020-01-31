use std::io::{self, Read};

use byteorder::ReadBytesExt;

use crate::{num::read_itf8, preservation_map::PreservationMap};

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
    let mut key_buf = vec![0; 2];

    for _ in 0..map_len {
        buf_reader.read_exact(&mut key_buf)?;

        match &key_buf[..] {
            b"RN" => {
                *map.read_names_included_mut() = read_bool(&mut buf_reader)?;
            }
            b"AP" => {
                *map.ap_data_series_delta_mut() = read_bool(&mut buf_reader)?;
            }
            b"RR" => {
                *map.read_names_included_mut() = read_bool(&mut buf_reader)?;
            }
            b"SM" => {
                let buf = map.substitution_matrix_mut();
                reader.read_exact(&mut buf[..])?;
            }
            b"TD" => {
                let td_len = read_itf8(reader)?;
                let buf = map.tag_ids_dictionary_mut();
                buf.resize(td_len as usize, Default::default());
                reader.read_exact(buf)?;
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid map key",
                ))
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
