mod reference_sequence_names;

use std::io::{self, Write};

use self::reference_sequence_names::write_reference_sequence_names;
use crate::{binning_index::index::Header, io::writer::num::write_i32_le};

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

    let col_end = header.end_position_index().map_or(Ok(0), |mut i| {
        i = i.checked_add(1).expect("attempt to add with overflow");
        i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;
    write_i32_le(writer, col_end)?;

    let meta = i32::from(header.line_comment_prefix());
    write_i32_le(writer, meta)?;

    let skip = i32::try_from(header.line_skip_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, skip)?;

    write_reference_sequence_names(writer, header.reference_sequence_names())?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;

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
}
