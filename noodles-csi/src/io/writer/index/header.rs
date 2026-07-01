mod reference_sequence_names;

use std::io::{self, Write};

use self::reference_sequence_names::write_reference_sequence_names;
use crate::{
    binning_index::index::{Header, header::Format},
    io::writer::num::write_i32_le,
};

pub(super) fn write_aux<W>(writer: &mut W, header: Option<&Header>) -> io::Result<()>
where
    W: Write,
{
    let mut aux = Vec::new();

    if let Some(hdr) = header {
        write_header(&mut aux, hdr)?;
    }

    let l_aux =
        i32::try_from(aux.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, l_aux)?;

    writer.write_all(&aux)?;

    Ok(())
}

pub(crate) fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: Write,
{
    let format = i32::from(header.format());
    write_i32_le(writer, format)?;

    let reference_sequence_name_index = header
        .reference_sequence_name_index()
        .checked_add(1)
        .expect("attempt to add with overflow");
    let col_seq = i32::try_from(reference_sequence_name_index)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, col_seq)?;

    let start_position_index = header
        .start_position_index()
        .checked_add(1)
        .expect("attempt to add with overflow");
    let col_beg = i32::try_from(start_position_index)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, col_beg)?;

    write_end_position_index(
        writer,
        header.format(),
        header.start_position_index(),
        header.end_position_index(),
    )?;

    let meta = i32::from(header.line_comment_prefix());
    write_i32_le(writer, meta)?;

    let skip = i32::try_from(header.line_skip_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, skip)?;

    write_reference_sequence_names(writer, header.reference_sequence_names())?;

    Ok(())
}

fn write_end_position_index<W>(
    writer: &mut W,
    format: Format,
    start_position_index: usize,
    end_position_index: Option<usize>,
) -> io::Result<()>
where
    W: Write,
{
    const SPECIALIZED_END_VALUE: i32 = 0;

    if matches!(format, Format::Sam | Format::Vcf) {
        if end_position_index.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid end position index for format",
            ));
        } else {
            write_i32_le(writer, SPECIALIZED_END_VALUE)?;
        }
    } else {
        let i = end_position_index.unwrap_or(start_position_index);
        let j = i.checked_add(1).expect("attempt to add with overflow");
        let n = i32::try_from(j).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_i32_le(writer, n)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;
    use crate::binning_index::index::header::format::CoordinateSystem;

    #[test]
    fn test_write_aux() -> io::Result<()> {
        let mut buf = Vec::new();
        write_aux(&mut buf, None)?;
        let expected = [0x00, 0x00, 0x00, 0x00];
        assert_eq!(buf, expected);

        let names = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();
        let header = crate::binning_index::index::header::Builder::vcf()
            .set_reference_sequence_names(names)
            .build();

        buf.clear();
        write_aux(&mut buf, Some(&header))?;

        let expected = [
            0x24, 0x00, 0x00, 0x00, // l_aux = 36
            0x02, 0x00, 0x00, 0x00, // format = 2 (VCF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1 (1-based)
            0x02, 0x00, 0x00, 0x00, // col_beg = 2 (1-based)
            0x00, 0x00, 0x00, 0x00, // col_end = None (1-based)
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
            0x08, 0x00, 0x00, 0x00, // l_nm = 8
            b's', b'q', b'0', 0x00, // names[0] = "sq0"
            b's', b'q', b'1', 0x00, // names[1] = "sq1"
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_end_position_index() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            format: Format,
            start_position_index: usize,
            end_position_index: Option<usize>,
            expected: i32,
        ) -> io::Result<()> {
            buf.clear();
            write_end_position_index(buf, format, start_position_index, end_position_index)?;
            assert_eq!(buf, &expected.to_le_bytes());
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Format::Sam, 5, None, 0)?;
        t(&mut buf, Format::Vcf, 5, None, 0)?;
        t(&mut buf, Format::Generic(CoordinateSystem::Gff), 5, None, 6)?;
        t(
            &mut buf,
            Format::Generic(CoordinateSystem::Gff),
            5,
            Some(8),
            9,
        )?;

        buf.clear();
        assert!(matches!(
            write_end_position_index(&mut buf, Format::Sam, 5, Some(8)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        assert!(matches!(
            write_end_position_index(&mut buf, Format::Vcf, 5, Some(8)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
