use std::{
    convert::TryFrom,
    io::{self, BufRead, Read},
};

use byteorder::ReadBytesExt;
use noodles_bam as bam;

use crate::{
    container::compression_header::{
        preservation_map::Key, PreservationMap, SubstitutionMatrix, TagIdsDictionary,
    },
    num::read_itf8,
    record,
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

    let mut read_names_included = true;
    let mut ap_data_series_delta = true;
    let mut reference_required = true;
    let mut substitution_matrix = None;
    let mut tag_ids_dictionary = None;

    let mut key_buf = [0; 2];

    for _ in 0..map_len {
        buf_reader.read_exact(&mut key_buf)?;

        let key = Key::try_from(&key_buf[..])
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

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
    SubstitutionMatrix::try_from(&buf[..])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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

        let mut line = Vec::new();

        for chunk in keys_buf.chunks_exact(3) {
            let (t0, t1, ty) = (chunk[0], chunk[1], chunk[2]);

            let tag = [t0, t1];
            let ty = bam::record::data::field::value::Type::try_from(ty)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let key = record::tag::Key::new(tag, ty);

            line.push(key);
        }

        dictionary.push(line);
    }

    Ok(dictionary)
}
