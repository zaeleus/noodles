use std::io::{self, Write};

use super::write_newline;
use crate::record::Sequence;

pub(super) fn write_sequence<W>(
    writer: &mut W,
    sequence: &Sequence,
    line_bases: usize,
) -> io::Result<()>
where
    W: Write,
{
    for bases in sequence.as_ref().chunks(line_bases) {
        writer.write_all(bases)?;
        write_newline(writer)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_sequence() -> io::Result<()> {
        let mut writer = Vec::new();
        let sequence = Sequence::from(b"AC".to_vec());
        write_sequence(&mut writer, &sequence, 4)?;
        assert_eq!(writer, b"AC\n");

        writer.clear();
        let sequence = Sequence::from(b"ACGT".to_vec());
        write_sequence(&mut writer, &sequence, 4)?;
        assert_eq!(writer, b"ACGT\n");

        writer.clear();
        let sequence = Sequence::from(b"ACGTACGT".to_vec());
        write_sequence(&mut writer, &sequence, 4)?;
        assert_eq!(writer, b"ACGT\nACGT\n");

        writer.clear();
        let sequence = Sequence::from(b"ACGTACGTAC".to_vec());
        write_sequence(&mut writer, &sequence, 4)?;
        assert_eq!(writer, b"ACGT\nACGT\nAC\n");

        Ok(())
    }
}
