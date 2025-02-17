use std::io::{self, Write};

use crate::{
    container::compression_header::preservation_map::TagSets,
    io::writer::{num::write_itf8, Record},
};

const NUL: u8 = 0x00;

pub(super) fn write_tag_sets<W>(writer: &mut W, tag_sets: &TagSets) -> io::Result<()>
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

pub(super) fn build_tag_sets(records: &[Record]) -> TagSets {
    use crate::container::compression_header::preservation_map::tag_sets::Key;

    let mut tag_sets = Vec::new();

    for record in records {
        let set = record
            .data
            .iter()
            .map(|(tag, value)| Key::new(*tag, value.ty()))
            .collect();

        if tag_sets.iter().any(|s| s == &set) {
            continue;
        } else {
            tag_sets.push(set);
        }
    }

    tag_sets
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::data::field::{Tag, Type};

    use super::*;
    use crate::container::compression_header::preservation_map::tag_sets::Key;

    #[test]
    fn test_write_tag_sets() -> io::Result<()> {
        let mut buf = Vec::new();

        let tag_sets = vec![
            vec![Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8)],
            vec![
                Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8),
                Key::new(Tag::COMMENT, Type::String),
            ],
        ];

        write_tag_sets(&mut buf, &tag_sets)?;

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
