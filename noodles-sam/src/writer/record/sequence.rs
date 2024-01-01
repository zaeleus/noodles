use std::io::{self, Write};

use super::MISSING;
use crate::alignment::record_buf::Sequence;

pub fn write_sequence<W>(writer: &mut W, read_length: usize, sequence: &Sequence) -> io::Result<()>
where
    W: Write,
{
    let sequence = sequence.as_ref();

    // § 1.4.10 "`SEQ`" (2022-08-22): "This field can be a '*' when the sequence is not stored."
    if sequence.is_empty() {
        writer.write_all(&[MISSING])?;
        return Ok(());
    }

    // § 1.4.10 "`SEQ`" (2022-08-22): "If not a '*', the length of the sequence must equal the sum
    // of lengths of `M`/`I`/`S`/`=`/`X` operations in `CIGAR`."
    if read_length > 0 && sequence.len() != read_length {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "read length-sequence length mismatch",
        ));
    }

    if !is_valid(sequence) {
        return Err(io::Error::from(io::ErrorKind::InvalidInput));
    }

    writer.write_all(sequence)?;

    Ok(())
}

fn is_valid(sequence: &[u8]) -> bool {
    sequence
        .iter()
        .all(|&b| matches!(b, b'A'..=b'Z' | b'a'..=b'z' | b'=' | b'.'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        write_sequence(&mut buf, 0, &Sequence::default())?;
        assert_eq!(buf, b"*");

        buf.clear();
        let sequence = Sequence::from(b"ACGT".to_vec());
        write_sequence(&mut buf, 4, &sequence)?;
        assert_eq!(buf, b"ACGT");

        buf.clear();
        let sequence = Sequence::from(b"ACGT".to_vec());
        write_sequence(&mut buf, 0, &sequence)?;
        assert_eq!(buf, b"ACGT");

        buf.clear();
        assert!(matches!(
            write_sequence(&mut buf, 1, &sequence),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        let sequence = Sequence::from(vec![b'!']);
        buf.clear();
        assert!(matches!(
            write_sequence(&mut buf, 1, &sequence),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }
}
