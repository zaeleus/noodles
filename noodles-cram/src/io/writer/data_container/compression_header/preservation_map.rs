use std::io::{self, Write};

use byteorder::WriteBytesExt;

use crate::{
    data_container::compression_header::{
        preservation_map::{Key, SubstitutionMatrix, TagSets},
        PreservationMap,
    },
    io::writer::{collections::write_array, num::write_itf8},
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
    let buf = encode(preservation_map)?;
    write_array(writer, &buf)
}

fn encode(preservation_map: &PreservationMap) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    encode_inner(&mut buf, preservation_map)?;
    Ok(buf)
}

fn encode_inner<W>(writer: &mut W, preservation_map: &PreservationMap) -> io::Result<()>
where
    W: Write,
{
    write_itf8(writer, MAP_LENGTH)?;

    write_key(writer, Key::ReadNamesIncluded)?;
    write_bool(writer, preservation_map.read_names_included())?;

    write_key(writer, Key::ApDataSeriesDelta)?;
    write_bool(writer, preservation_map.ap_data_series_delta())?;

    write_key(writer, Key::ReferenceRequired)?;
    write_bool(writer, preservation_map.is_reference_required())?;

    write_key(writer, Key::SubstitutionMatrix)?;
    write_substitution_matrix(writer, preservation_map.substitution_matrix())?;

    write_key(writer, Key::TagSets)?;
    write_tag_sets(writer, preservation_map.tag_sets())?;

    Ok(())
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

fn write_tag_sets<W>(writer: &mut W, tag_sets: &TagSets) -> io::Result<()>
where
    W: Write,
{
    use noodles_bam::record::codec::encoder::data::field::ty::encode;

    let mut buf = Vec::new();

    for keys in tag_sets.iter() {
        for key in keys {
            let tag = key.tag();
            buf.extend_from_slice(tag.as_ref());

            let ty = key.ty();
            buf.push(encode(ty));
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
    use noodles_sam::alignment::record::data::field::{Tag, Type};

    use super::*;
    use crate::data_container::compression_header::preservation_map::tag_sets::Key;

    #[test]
    fn test_write_tag_ids_dictionary() -> io::Result<()> {
        let mut buf = Vec::new();

        let tag_ids_dictionary = TagSets::from(vec![
            vec![Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8)],
            vec![
                Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8),
                Key::new(Tag::COMMENT, Type::String),
            ],
        ]);

        write_tag_sets(&mut buf, &tag_ids_dictionary)?;

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
