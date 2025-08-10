use std::io::{self, Write};

use bstr::{BStr, ByteSlice};
use noodles_csi::binning_index::index::header::ReferenceSequenceNames;

use crate::io::writer::num::write_i32_le;

const NUL: u8 = 0x00;

pub(super) fn write_reference_sequence_names<W>(
    writer: &mut W,
    reference_sequence_names: &ReferenceSequenceNames,
) -> io::Result<()>
where
    W: Write,
{
    // Add 1 for each trailing NUL.
    let len = reference_sequence_names
        .iter()
        .map(|n| n.len() + 1)
        .sum::<usize>();
    let l_nm = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, l_nm)?;

    for reference_sequence_name in reference_sequence_names {
        if !is_valid(reference_sequence_name.as_ref()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid reference sequence name",
            ));
        }

        writer.write_all(reference_sequence_name)?;
        writer.write_all(&[NUL])?;
    }

    Ok(())
}

fn is_valid(s: &BStr) -> bool {
    s.find_byte(NUL).is_none()
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;

    #[test]
    fn test_write_reference_sequence_names() -> io::Result<()> {
        let mut buf = Vec::new();

        let reference_sequence_names = ReferenceSequenceNames::default();
        buf.clear();
        write_reference_sequence_names(&mut buf, &reference_sequence_names)?;
        assert_eq!(buf, [0x00, 0x00, 0x00, 0x00]); // l_nm = 0

        let reference_sequence_names = [BString::from("sq0")].into_iter().collect();
        buf.clear();
        write_reference_sequence_names(&mut buf, &reference_sequence_names)?;
        assert_eq!(
            buf,
            [
                0x04, 0x00, 0x00, 0x00, // l_nm = 4
                b's', b'q', b'0', 0x00, // names[0] = "sq0"
            ]
        );

        let reference_sequence_names = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();
        buf.clear();
        write_reference_sequence_names(&mut buf, &reference_sequence_names)?;
        assert_eq!(
            buf,
            [
                0x08, 0x00, 0x00, 0x00, // l_nm = 8
                b's', b'q', b'0', 0x00, // names[0] = "sq0"
                b's', b'q', b'1', 0x00, // names[1] = "sq1"
            ]
        );

        let reference_sequence_names = [BString::from("sq\x000")].into_iter().collect();
        buf.clear();
        assert!(matches!(
            write_reference_sequence_names(&mut buf, &reference_sequence_names),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
