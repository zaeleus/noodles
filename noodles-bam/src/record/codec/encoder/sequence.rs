use std::io;

use bytes::BufMut;
use noodles_sam::alignment::record_buf::Sequence;

pub fn put_sequence<B>(dst: &mut B, read_length: usize, sequence: &Sequence) -> io::Result<()>
where
    B: BufMut,
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

    let mut bases = sequence.as_ref().iter().copied();

    while let Some(l) = bases.next() {
        // § 4.2.3 "SEQ and QUAL encoding" (2021-06-03): "When `l_seq` is odd the bottom 4 bits of
        // the last byte are undefined, but we recommend writing these as zero."
        let r = bases.next().unwrap_or(EQ);
        let b = encode_base(l) << 4 | encode_base(r);
        dst.put_u8(b);
    }

    Ok(())
}

fn encode_base(n: u8) -> u8 {
    match n {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        // § 4.2.3 SEQ and QUAL encoding (2021-06-03): "The case-insensitive base codes ... are
        // mapped to [0, 15] respectively with all other characters mapping to 'N' (value 15)".
        _ => 15,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_sequence() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, sequence: &Sequence, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_sequence(buf, sequence.len(), sequence)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Sequence::default(), &[])?;
        t(&mut buf, &Sequence::from(b"ACG".to_vec()), &[0x12, 0x40])?;
        t(&mut buf, &Sequence::from(b"ACGT".to_vec()), &[0x12, 0x48])?;

        buf.clear();
        put_sequence(&mut buf, 2, &Sequence::default())?;
        assert!(buf.is_empty());

        buf.clear();
        let sequence = Sequence::from(b"A".to_vec());
        assert!(matches!(
            put_sequence(&mut buf, 2, &sequence),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'='), 0);
        assert_eq!(encode_base(b'A'), 1);
        assert_eq!(encode_base(b'C'), 2);
        assert_eq!(encode_base(b'M'), 3);
        assert_eq!(encode_base(b'G'), 4);
        assert_eq!(encode_base(b'R'), 5);
        assert_eq!(encode_base(b'S'), 6);
        assert_eq!(encode_base(b'V'), 7);
        assert_eq!(encode_base(b'T'), 8);
        assert_eq!(encode_base(b'W'), 9);
        assert_eq!(encode_base(b'Y'), 10);
        assert_eq!(encode_base(b'H'), 11);
        assert_eq!(encode_base(b'K'), 12);
        assert_eq!(encode_base(b'D'), 13);
        assert_eq!(encode_base(b'B'), 14);
        assert_eq!(encode_base(b'N'), 15);

        assert_eq!(encode_base(b'X'), 15);
    }
}
