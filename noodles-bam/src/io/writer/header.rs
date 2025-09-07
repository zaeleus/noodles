use std::{
    ffi::CString,
    io::{self, Write},
    num::NonZero,
};

use noodles_sam::{self as sam, header::ReferenceSequences};

use crate::io::writer::num::{write_i32_le, write_u32_le};

pub(super) fn write_header<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: Write,
{
    write_raw_header(writer, header)?;
    write_reference_sequences(writer, header.reference_sequences())?;
    Ok(())
}

fn write_raw_header<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: Write,
{
    use crate::io::MAGIC_NUMBER;

    writer.write_all(&MAGIC_NUMBER)?;

    let text = serialize_header(header)?;
    let l_text =
        i32::try_from(text.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, l_text)?;

    writer.write_all(&text)?;

    Ok(())
}

fn serialize_header(header: &sam::Header) -> io::Result<Vec<u8>> {
    let mut writer = sam::io::Writer::new(Vec::new());
    writer.write_header(header)?;
    Ok(writer.into_inner())
}

pub fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &ReferenceSequences,
) -> io::Result<()>
where
    W: Write,
{
    let n_ref = i32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, n_ref)?;

    for (name, reference_sequence) in reference_sequences {
        write_reference_sequence(writer, name, reference_sequence.length())?;
    }

    Ok(())
}

fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence_name: &[u8],
    length: NonZero<usize>,
) -> io::Result<()>
where
    W: Write,
{
    let c_name = CString::new(reference_sequence_name)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let name = c_name.as_bytes_with_nul();

    let l_name =
        u32::try_from(name.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, l_name)?;
    writer.write_all(name)?;

    let l_ref = i32::try_from(usize::from(length))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, l_ref)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;

    #[test]
    fn test_write_raw_header() -> Result<(), Box<dyn std::error::Error>> {
        use sam::header::record::value::{
            Map,
            map::{self, header::Version},
        };

        let header = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        let mut buf = Vec::new();
        write_raw_header(&mut buf, &header)?;

        let mut expected = vec![
            b'B', b'A', b'M', 0x01, // magic
            0x0b, 0x00, 0x00, 0x00, // l_text = 11
        ];
        expected.extend_from_slice(b"@HD\tVN:1.6\n"); // text

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_reference_sequences() -> io::Result<()> {
        use sam::header::record::value::{Map, map::ReferenceSequence};

        let reference_sequences = [(
            BString::from("sq0"),
            Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
        )]
        .into_iter()
        .collect();

        let mut buf = Vec::new();
        write_reference_sequences(&mut buf, &reference_sequences)?;

        let expected = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            b's', b'q', b'0', 0x00, // ref[0].name = b"sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_reference_sequence() -> io::Result<()> {
        let mut buf = Vec::new();
        write_reference_sequence(&mut buf, b"sq0", const { NonZero::new(8).unwrap() })?;

        let expected = [
            0x04, 0x00, 0x00, 0x00, // l_name = 4
            0x73, 0x71, 0x30, 0x00, // name = b"sq0\x00"
            0x08, 0x00, 0x00, 0x00, // l_ref = 8
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
