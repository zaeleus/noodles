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
    let src = sequence.as_ref();
    let base_count = sequence.len();

    let pair_count = base_count / 2;
    let (pairs, rest) = src.split_at(src.len().min(pair_count));

    for &n in pairs {
        let bases = decode_bases(n);
        writer.write_all(&bases)?;
    }

    if !base_count.is_multiple_of(2)
        && let Some(&n) = rest.first()
    {
        let [b, _] = decode_bases(n);
        writer.write_all(&[b])?;
    }

    Ok(())
}

const CODES: [[u8; 2]; 256] = build_codes();

const fn build_codes() -> [[u8; 2]; 256] {
    const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

    let mut table = [[0u8; 2]; 256];
    let mut i = 0;

    while i < 256 {
        table[i] = [BASES[i >> 4], BASES[i & 0xf]];
        i += 1;
    }

    table
}

fn decode_bases(n: u8) -> [u8; 2] {
    CODES[usize::from(n)]
}

fn write_raw_sequence<W>(writer: &mut W, sequence: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if sequence.iter().all(|&b| is_valid_base(b)) {
        writer.write_all(sequence)
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidInput))
    }
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
        fn t(buf: &mut Vec<u8>, sequence: &FourBitPacked, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_four_bit_packed_sequence(buf, sequence)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &FourBitPacked::new(&[0x12, 0x40], 3), b"ACG")?;
        t(&mut buf, &FourBitPacked::new(&[0x12, 0x48], 4), b"ACGT")?;
        t(
            &mut buf,
            &FourBitPacked::new(&[0x12, 0x48, 0x00], 4),
            b"ACGT",
        )?;

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
