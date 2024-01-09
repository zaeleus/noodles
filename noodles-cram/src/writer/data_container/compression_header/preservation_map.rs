use std::io::{self, Write};

use byteorder::WriteBytesExt;

use crate::{
    data_container::compression_header::{
        preservation_map::{Key, SubstitutionMatrix, TagIdsDictionary},
        PreservationMap,
    },
    writer::num::write_itf8,
};

const MAP_LENGTH: i32 = 5;

const FALSE: u8 = 0x00;
const TRUE: u8 = 0x01;

const NUL: u8 = 0x00;

pub(super) fn write_preservation_map<W>(
    writer: &mut W,
    preservation_map: &PreservationMap,
) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    write_itf8(&mut buf, MAP_LENGTH)?;

    write_key(&mut buf, Key::ReadNamesIncluded)?;
    write_bool(&mut buf, preservation_map.read_names_included())?;

    write_key(&mut buf, Key::ApDataSeriesDelta)?;
    write_bool(&mut buf, preservation_map.ap_data_series_delta())?;

    write_key(&mut buf, Key::ReferenceRequired)?;
    write_bool(&mut buf, preservation_map.is_reference_required())?;

    write_key(&mut buf, Key::SubstitutionMatrix)?;
    write_substitution_matrix(&mut buf, preservation_map.substitution_matrix())?;

    write_key(&mut buf, Key::TagIdsDictionary)?;
    write_tag_ids_dictionary(&mut buf, preservation_map.tag_ids_dictionary())?;

    let data_len =
        i32::try_from(buf.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, data_len)?;

    writer.write_all(&buf)
}

fn write_key<W>(writer: &mut W, key: Key) -> io::Result<()>
where
    W: Write,
{
    let data = <[u8; 2]>::from(key);
    writer.write_all(&data)
}

fn write_bool<W>(writer: &mut W, value: bool) -> io::Result<()>
where
    W: Write,
{
    if value {
        writer.write_u8(TRUE)
    } else {
        writer.write_u8(FALSE)
    }
}

fn write_substitution_matrix<W>(
    writer: &mut W,
    substitution_matrix: &SubstitutionMatrix,
) -> io::Result<()>
where
    W: Write,
{
    let buf = <[u8; 5]>::from(substitution_matrix);
    writer.write_all(&buf)
}

fn write_tag_ids_dictionary<W>(
    writer: &mut W,
    tag_ids_dictionary: &TagIdsDictionary,
) -> io::Result<()>
where
    W: Write,
{
    use noodles_bam::record::codec::encoder::data::field::ty::type_to_u8;

    let mut buf = Vec::new();

    for keys in tag_ids_dictionary.iter() {
        for key in keys {
            let tag = key.tag();
            buf.extend_from_slice(tag.as_ref());

            let ty = key.ty();
            buf.push(type_to_u8(ty));
        }

        buf.push(NUL);
    }

    let data_len =
        i32::try_from(buf.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, data_len)?;
    writer.write_all(&buf)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_sam::{alignment::record::data::field::Type, record::data::field::tag};

    use super::*;
    use crate::data_container::compression_header::preservation_map::tag_ids_dictionary::Key;

    #[test]
    fn test_write_tag_ids_dictionary() -> io::Result<()> {
        let mut buf = Vec::new();

        let tag_ids_dictionary = TagIdsDictionary::from(vec![
            vec![Key::new(tag::ALIGNMENT_HIT_COUNT, Type::Int8)],
            vec![
                Key::new(tag::ALIGNMENT_HIT_COUNT, Type::Int8),
                Key::new(tag::COMMENT, Type::String),
            ],
        ]);

        write_tag_ids_dictionary(&mut buf, &tag_ids_dictionary)?;

        let expected = [
            0x0b, // data_len
            0x4e, 0x48, 0x63, // NH:c
            0x00, // nul
            0x4e, 0x48, 0x63, // NH:c
            0x43, 0x4f, 0x5a, // CO:Z
            0x00, // nul
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
