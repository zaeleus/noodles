mod reference_sequence_names;

use std::io::{self, Write};

use noodles_csi::binning_index::index::Header;

use self::reference_sequence_names::write_reference_sequence_names;
use crate::io::writer::num::write_i32_le;

pub(super) fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
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

    write_end_position_index(writer, header.end_position_index())?;

    let meta = i32::from(header.line_comment_prefix());
    write_i32_le(writer, meta)?;

    let skip = i32::try_from(header.line_skip_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, skip)?;

    write_reference_sequence_names(writer, header.reference_sequence_names())?;

    Ok(())
}

fn write_end_position_index<W>(writer: &mut W, i: Option<usize>) -> io::Result<()>
where
    W: Write,
{
    let n = i.map_or(Ok(0), |mut j| {
        j = j.checked_add(1).expect("attempt to add with overflow");
        i32::try_from(j).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;

    write_i32_le(writer, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_end_position_index() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_end_position_index(&mut buf, None)?;
        assert_eq!(buf, 0i32.to_le_bytes());

        buf.clear();
        write_end_position_index(&mut buf, Some(0))?;
        assert_eq!(buf, 1i32.to_le_bytes());

        Ok(())
    }
}
