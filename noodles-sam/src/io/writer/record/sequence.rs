use std::io::{self, Write};

use super::MISSING;
use crate::alignment::record::{Sequence, SequenceRef, sequence_ref::FourBitPacked};

pub(super) fn write_sequence<W>(
    writer: &mut W,
    read_length: usize,
    sequence: SequenceRef<'_>,
) -> io::Result<()>
where
    W: Write,
{
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

    match sequence {
        SequenceRef::FourBitPacked(sequence) => write_four_bit_packed_sequence(writer, &sequence)?,
        SequenceRef::Raw(sequence) => write_raw_sequence(writer, sequence)?,
        SequenceRef::Sequence(sequence) => write_generic_sequence(writer, sequence)?,
    }

    Ok(())
}

fn write_four_bit_packed_sequence<W>(writer: &mut W, sequence: &FourBitPacked) -> io::Result<()>
where
    W: Write,
{
    for base in sequence.iter() {
        // SAFETY: `base` is guaranteed to be a valid base.
        writer.write_all(&[base])?;
    }

    Ok(())
}

fn write_raw_sequence<W>(writer: &mut W, sequence: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if !sequence.iter().all(|&b| is_valid_base(b)) {
        return Err(io::Error::from(io::ErrorKind::InvalidInput));
    }

    for &b in sequence {
        writer.write_all(&[b])?;
    }

    Ok(())
}

fn write_generic_sequence<W, S>(writer: &mut W, sequence: S) -> io::Result<()>
where
    W: Write,
    S: Sequence,
{
    for base in sequence.iter() {
        if !is_valid_base(base) {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        writer.write_all(&[base])?;
    }

    Ok(())
}

// § 1.4 "The alignment section: mandatory fields" (2024-11-06): `[A-Za-z=.]+`.
fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A'..=b'Z' | b'a'..=b'z' | b'=' | b'.')
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::record_buf::Sequence as SequenceBuf;

    #[test]
    fn test_write_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let sequence = SequenceBuf::default();
        let s = SequenceRef::Sequence(Box::new(&sequence));
        write_sequence(&mut buf, 0, s)?;
        assert_eq!(buf, b"*");

        buf.clear();
        let sequence = SequenceBuf::from(b"ACGT");
        let s = SequenceRef::Sequence(Box::new(&sequence));
        write_sequence(&mut buf, 4, s)?;
        assert_eq!(buf, b"ACGT");

        buf.clear();
        let sequence = SequenceBuf::from(b"ACGT");
        let s = SequenceRef::Sequence(Box::new(&sequence));
        write_sequence(&mut buf, 0, s)?;
        assert_eq!(buf, b"ACGT");

        buf.clear();
        let sequence = SequenceBuf::from(b"ACGT");
        let s = SequenceRef::Sequence(Box::new(&sequence));
        assert!(matches!(
            write_sequence(&mut buf, 1, s),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        buf.clear();
        let sequence = SequenceBuf::from(vec![b'!']);
        let s = SequenceRef::Sequence(Box::new(&sequence));
        assert!(matches!(
            write_sequence(&mut buf, 1, s),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_write_four_bit_packed_sequence() -> io::Result<()> {
        let mut buf = Vec::new();
        let sequence = FourBitPacked::new(&[0x12, 0x48], 4); // ACGT
        write_four_bit_packed_sequence(&mut buf, &sequence)?;
        assert_eq!(buf, b"ACGT");
        Ok(())
    }

    #[test]
    fn test_write_raw_sequence() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_raw_sequence(&mut buf, b"ACGT")?;
        assert_eq!(buf, b"ACGT");

        buf.clear();
        assert!(matches!(
            write_raw_sequence(&mut buf, &[0xf0, 0x9f, 0x8d, 0x9c]),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid_base() {
        for b in (b'A'..=b'Z').chain(b'a'..=b'z').chain(*b"=.") {
            assert!(is_valid_base(b));
        }

        assert!(!is_valid_base(b'!'));
    }
}
