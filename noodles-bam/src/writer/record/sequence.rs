use bytes::BufMut;
use noodles_sam::record::{sequence::Base, Sequence};

pub fn put_sequence<B>(dst: &mut B, sequence: &Sequence)
where
    B: BufMut,
{
    let mut bases = sequence.as_ref().iter().copied();

    while let Some(l) = bases.next() {
        // § 4.2.3 "SEQ and QUAL encoding" (2021-06-03): "When `l_seq` is odd the bottom 4 bits of
        // the last byte are undefined, but we recommend writing these as zero."
        let r = bases.next().unwrap_or(Base::Eq);
        let b = encode_base(l) << 4 | encode_base(r);
        dst.put_u8(b);
    }
}

fn encode_base(base: Base) -> u8 {
    match base {
        Base::Eq => 0,
        Base::A => 1,
        Base::C => 2,
        Base::M => 3,
        Base::G => 4,
        Base::R => 5,
        Base::S => 6,
        Base::V => 7,
        Base::T => 8,
        Base::W => 9,
        Base::Y => 10,
        Base::H => 11,
        Base::K => 12,
        Base::D => 13,
        Base::B => 14,
        // § 4.2.3 SEQ and QUAL encoding (2021-06-03): "The case-insensitive base codes ... are
        // mapped to [0, 15] respectively with all other characters mapping to 'N' (value 15)".
        _ => 15,
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam as sam;

    use super::*;

    #[test]
    fn test_put_sequence() -> Result<(), sam::record::sequence::ParseError> {
        use sam::record::Sequence;

        fn t(buf: &mut Vec<u8>, sequence: &Sequence, expected: &[u8]) {
            buf.clear();
            put_sequence(buf, sequence);
            assert_eq!(buf, expected);
        }

        let mut buf = Vec::new();

        t(&mut buf, &Sequence::default(), &[]);
        t(&mut buf, &"ACG".parse()?, &[0x12, 0x40]);
        t(&mut buf, &"ACGT".parse()?, &[0x12, 0x48]);

        Ok(())
    }

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(Base::Eq), 0);
        assert_eq!(encode_base(Base::A), 1);
        assert_eq!(encode_base(Base::C), 2);
        assert_eq!(encode_base(Base::M), 3);
        assert_eq!(encode_base(Base::G), 4);
        assert_eq!(encode_base(Base::R), 5);
        assert_eq!(encode_base(Base::S), 6);
        assert_eq!(encode_base(Base::V), 7);
        assert_eq!(encode_base(Base::T), 8);
        assert_eq!(encode_base(Base::W), 9);
        assert_eq!(encode_base(Base::Y), 10);
        assert_eq!(encode_base(Base::H), 11);
        assert_eq!(encode_base(Base::K), 12);
        assert_eq!(encode_base(Base::D), 13);
        assert_eq!(encode_base(Base::B), 14);
        assert_eq!(encode_base(Base::N), 15);

        assert_eq!(encode_base(Base::X), 15);
    }
}
