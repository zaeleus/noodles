mod reference_sequence_names;

use std::io::{self, Write};

use noodles_csi::binning_index::index::{Header, header::Format};

use self::reference_sequence_names::write_reference_sequence_names;
use crate::io::writer::num::write_i32_le;

pub(super) fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
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
    use noodles_csi::binning_index::index::header::format::CoordinateSystem;

    use super::*;

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
