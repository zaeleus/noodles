use std::io::{self, Write};

use crate::{
    container::compression_header::preservation_map::{SubstitutionMatrix, TagIdsDictionary},
    num::{write_itf8, Itf8},
};

const NUL: u8 = 0x00;

#[allow(dead_code)]
fn write_substituion_matrix<W>(
    writer: &mut W,
    substitution_matrix: &SubstitutionMatrix,
) -> io::Result<()>
where
    W: Write,
{
    // FIXME: unnecessary clone
    let buf = <[u8; 5]>::from(substitution_matrix.clone());
    writer.write_all(&buf)
}

#[allow(dead_code)]
fn write_tag_ids_dictionary<W>(
    writer: &mut W,
    tag_ids_dictionary: &TagIdsDictionary,
) -> io::Result<()>
where
    W: Write,
{
    let mut buf = Vec::new();

    for keys in tag_ids_dictionary.iter() {
        for key in keys {
            let id = key.id();
            buf.push((id >> 16) as u8);
            buf.push((id >> 8) as u8);
            buf.push(id as u8);
        }

        buf.push(NUL);
    }

    // FIXME: usize => Itf8 cast
    let data_len = buf.len() as Itf8;
    write_itf8(writer, data_len)?;
    writer.write_all(&buf)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_bam::record::data::field::value::Type;

    use crate::record::tag::Key;

    use super::*;

    #[test]
    fn test_write_tag_ids_dictionary() -> io::Result<()> {
        let mut buf = Vec::new();

        let tag_ids_dictionary = TagIdsDictionary::from(vec![
            vec![Key::new([b'N', b'H'], Type::Int8)],
            vec![
                Key::new([b'N', b'H'], Type::Int8),
                Key::new([b'C', b'O'], Type::String),
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
