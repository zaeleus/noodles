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

#[doc(hidden)]
pub fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: Write,
{
    write_format(writer, header.format())?;
    write_reference_sequence_name_index(writer, header.reference_sequence_name_index())?;
    write_start_position_index(writer, header.start_position_index())?;

    write_end_position_index(
        writer,
        header.format(),
        header.start_position_index(),
        header.end_position_index(),
    )?;

    write_line_comment_prefix(writer, header.line_comment_prefix())?;
    write_line_skip_count(writer, header.line_skip_count())?;

    write_reference_sequence_names(writer, header.reference_sequence_names())?;

    Ok(())
}

fn write_format<W>(writer: &mut W, format: Format) -> io::Result<()>
where
    W: Write,
{
    let n = i32::from(format);
    write_i32_le(writer, n)
}

fn write_reference_sequence_name_index<W>(writer: &mut W, i: usize) -> io::Result<()>
where
    W: Write,
{
    let j = i.checked_add(1).expect("attempt to add with overflow");
    let n = i32::try_from(j).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, n)?;
    Ok(())
}

fn write_start_position_index<W>(writer: &mut W, i: usize) -> io::Result<()>
where
    W: Write,
{
    let j = i.checked_add(1).expect("attempt to add with overflow");
    let n = i32::try_from(j).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, n)?;
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

fn write_line_comment_prefix<W>(writer: &mut W, line_comment_prefix: u8) -> io::Result<()>
where
    W: Write,
{
    let n = i32::from(line_comment_prefix);
    write_i32_le(writer, n)
}

fn write_line_skip_count<W>(writer: &mut W, line_skip_count: u32) -> io::Result<()>
where
    W: Write,
{
    let n = i32::try_from(line_skip_count)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_i32_le(writer, n)?;

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
    fn test_write_format() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, format: Format, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_format(buf, format)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            Format::Generic(CoordinateSystem::Gff),
            &[0x00, 0x00, 0x00, 0x00],
        )?;
        t(
            &mut buf,
            Format::Generic(CoordinateSystem::Bed),
            &[0x00, 0x00, 0x01, 0x00],
        )?;
        t(&mut buf, Format::Sam, &[0x01, 0x00, 0x00, 0x00])?;
        t(&mut buf, Format::Vcf, &[0x02, 0x00, 0x00, 0x00])?;

        Ok(())
    }

    #[test]
    fn test_write_reference_sequence_name_index() -> io::Result<()> {
        let mut buf = Vec::new();
        write_reference_sequence_name_index(&mut buf, 0)?;
        assert_eq!(buf, 1i32.to_le_bytes());
        Ok(())
    }

    #[test]
    fn test_write_start_position_index() -> io::Result<()> {
        let mut buf = Vec::new();
        write_start_position_index(&mut buf, 0)?;
        assert_eq!(buf, 1i32.to_le_bytes());
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

    #[test]
    fn test_write_line_comment_prefix() -> io::Result<()> {
        let mut buf = Vec::new();
        write_line_comment_prefix(&mut buf, b'#')?;
        assert_eq!(buf, i32::from(b'#').to_le_bytes());
        Ok(())
    }

    #[test]
    fn test_write_line_skip_count() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_line_skip_count(&mut buf, 8)?;
        assert_eq!(buf, 8i32.to_le_bytes());

        buf.clear();
        assert!(matches!(
            write_line_skip_count(&mut buf, u32::MAX),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
