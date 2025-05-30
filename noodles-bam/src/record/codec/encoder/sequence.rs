use std::io;

use noodles_sam::alignment::record::Sequence;

use super::num::{write_u8, write_u32_le};

pub(super) fn write_length(dst: &mut Vec<u8>, base_count: usize) -> io::Result<()> {
    let n =
        u32::try_from(base_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    write_u32_le(dst, n);

    Ok(())
}

pub(super) fn write_sequence<S>(
    dst: &mut Vec<u8>,
    read_length: usize,
    sequence: S,
) -> io::Result<()>
where
    S: Sequence,
{
    const EQ: u8 = b'=';

    if sequence.is_empty() {
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

    let mut bases = sequence.iter();

    while let Some(l) = bases.next() {
        // § 4.2.3 "SEQ and QUAL encoding" (2021-06-03): "When `l_seq` is odd the bottom 4 bits of
        // the last byte are undefined, but we recommend writing these as zero."
        let r = bases.next().unwrap_or(EQ);
        let n = (encode_base(l) << 4) | encode_base(r);
        write_u8(dst, n);
    }

    Ok(())
}

// § 4.2.3 "SEQ and QUAL encoding" (2023-11-16): "The case-insensitive base codes [...] are mapped
// to [0, 15] respectively with all other characters mapping to 'N' (value 15)".
fn encode_base(n: u8) -> u8 {
    const CODES: [u8; 256] = build_codes();
    CODES[usize::from(n)]
}

const fn build_codes() -> [u8; 256] {
    // § 4.2.3 "SEQ and QUAL encoding" (2024-11-06)
    const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
    const N: u8 = 0x0f;

    let mut i = 0;
    let mut codes = [N; 256];

    while i < BASES.len() {
        let base = BASES[i];

        // SAFETY: `i < BASES.len() <= u8::MAX`.
        let code = i as u8;

        // SAFETY: `base <= codes.len() <= u8::MAX`.
        codes[base as usize] = code;
        codes[base.to_ascii_lowercase() as usize] = code;

        i += 1;
    }

    codes
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record_buf::Sequence as SequenceBuf;

    use super::*;

    #[test]
    fn test_write_length() -> io::Result<()> {
        let mut buf = Vec::new();
        write_length(&mut buf, 8)?;
        assert_eq!(buf, [0x08, 0x00, 0x00, 0x00]);
        Ok(())
    }

    #[cfg(not(target_pointer_width = "16"))]
    #[test]
    fn test_write_length_with_out_of_range_base_count() {
        let mut buf = Vec::new();

        assert!(matches!(
            write_length(&mut buf, usize::MAX),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));
    }

    #[test]
    fn test_write_sequence() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, sequence: &SequenceBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_sequence(buf, sequence.len(), sequence)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &SequenceBuf::default(), &[])?;
        t(&mut buf, &SequenceBuf::from(b"ACG"), &[0x12, 0x40])?;
        t(&mut buf, &SequenceBuf::from(b"ACGT"), &[0x12, 0x48])?;

        buf.clear();
        write_sequence(&mut buf, 2, &SequenceBuf::default())?;
        assert!(buf.is_empty());

        buf.clear();
        let sequence = SequenceBuf::from(b"A");
        assert!(matches!(
            write_sequence(&mut buf, 2, &sequence),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_encode_base() {
        const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

        for (i, b) in (0..).zip(BASES) {
            assert_eq!(encode_base(b), i);
            assert_eq!(encode_base(b.to_ascii_lowercase()), i);
        }

        assert_eq!(encode_base(b'X'), 15);
        assert_eq!(encode_base(b'x'), 15);
    }
}
