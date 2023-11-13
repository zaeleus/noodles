use std::{
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::index::{header::ReferenceSequenceNames, Header};

pub(super) fn read_aux<R>(reader: &mut R) -> io::Result<Option<Header>>
where
    R: Read,
{
    let l_aux = reader.read_i32::<LittleEndian>().and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if l_aux > 0 {
        let mut aux_reader = reader.take(l_aux);
        read_header(&mut aux_reader).map(Some)
    } else {
        Ok(None)
    }
}

pub(crate) fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: Read,
{
    use crate::index::header::Format;

    let format = reader.read_i32::<LittleEndian>().and_then(|n| {
        Format::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let col_seq = reader.read_i32::<LittleEndian>().and_then(|i| {
        usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|n| {
                n.checked_sub(1)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid col_seq"))
            })
    })?;

    let col_beg = reader.read_i32::<LittleEndian>().and_then(|i| {
        usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|n| {
                n.checked_sub(1)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid col_beg"))
            })
    })?;

    let col_end = reader.read_i32::<LittleEndian>().and_then(|i| match i {
        0 => Ok(None),
        _ => usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|n| {
                n.checked_sub(1)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid col_end"))
            })
            .map(Some),
    })?;

    let meta = reader
        .read_i32::<LittleEndian>()
        .and_then(|b| u8::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let skip = reader.read_i32::<LittleEndian>().and_then(|n| {
        u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let names = read_names(reader)?;

    Ok(Header::builder()
        .set_format(format)
        .set_reference_sequence_name_index(col_seq)
        .set_start_position_index(col_beg)
        .set_end_position_index(col_end)
        .set_line_comment_prefix(meta)
        .set_line_skip_count(skip)
        .set_reference_sequence_names(names)
        .build())
}

fn read_names<R>(reader: &mut R) -> io::Result<ReferenceSequenceNames>
where
    R: Read,
{
    let l_nm = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut names = vec![0; l_nm];
    reader.read_exact(&mut names)?;

    parse_names(&names)
}

pub(crate) fn parse_names(mut src: &[u8]) -> io::Result<ReferenceSequenceNames> {
    const NUL: u8 = 0x00;

    let mut names = ReferenceSequenceNames::new();

    while let Some(i) = src.iter().position(|&b| b == NUL) {
        let (raw_name, rest) = src.split_at(i);

        let name =
            str::from_utf8(raw_name).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        if !names.insert(name.into()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate reference sequence name: {name}"),
            ));
        }

        src = &rest[1..];
    }

    if src.is_empty() {
        Ok(names)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid reference sequence names",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_aux() -> io::Result<()> {
        let src = [0x00, 0x00, 0x00, 0x00];
        let mut reader = &src[..];
        assert!(read_aux(&mut reader)?.is_none());

        let src = [
            0x20, 0x00, 0x00, 0x00, // l_aux = 32
            0x02, 0x00, 0x00, 0x00, // format = 2 (VCF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1 (1-based)
            0x02, 0x00, 0x00, 0x00, // col_beg = 2 (1-based)
            0x00, 0x00, 0x00, 0x00, // col_end = None (1-based)
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
            0x04, 0x00, 0x00, 0x00, // l_nm = 4
            b's', b'q', b'0', 0x00, // names = ["sq0"]
        ];
        let mut reader = &src[..];

        let actual = read_aux(&mut reader)?;
        let expected = crate::index::header::Builder::vcf()
            .set_reference_sequence_names([String::from("sq0")].into_iter().collect())
            .build();

        assert_eq!(actual, Some(expected));

        Ok(())
    }

    #[test]
    fn test_parse_names() -> io::Result<()> {
        let data = b"sq0\x00sq1\x00";
        let actual = parse_names(&data[..])?;
        let expected: ReferenceSequenceNames = [String::from("sq0"), String::from("sq1")]
            .into_iter()
            .collect();
        assert_eq!(actual, expected);

        let data = b"";
        assert!(parse_names(&data[..])?.is_empty());

        Ok(())
    }

    #[test]
    fn test_parse_names_with_duplicate_name() {
        let data = b"sq0\x00sq0\x00";

        assert!(matches!(
            parse_names(data),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }

    #[test]
    fn test_parse_names_with_trailing_data() {
        let data = b"sq0\x00sq1\x00sq2";

        assert!(matches!(
            parse_names(data),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
