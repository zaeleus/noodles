use std::io;

use bytes::Buf;
use noodles_sam::record::{sequence::Base, Sequence};

pub fn get_sequence<B>(src: &mut B, sequence: &mut Sequence, l_seq: usize) -> io::Result<()>
where
    B: Buf,
{
    let seq_len = (l_seq + 1) / 2;

    if src.remaining() < seq_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let seq = src.take(seq_len);
    let bases = seq
        .chunk()
        .iter()
        .flat_map(|&b| [decode_base(b >> 4), decode_base(b)]);

    sequence.clear();
    sequence.as_mut().extend(bases);
    sequence.as_mut().truncate(l_seq);

    src.advance(seq_len);

    Ok(())
}

pub fn decode_base(n: u8) -> Base {
    match n & 0x0f {
        0 => Base::Eq,
        1 => Base::A,
        2 => Base::C,
        3 => Base::M,
        4 => Base::G,
        5 => Base::R,
        6 => Base::S,
        7 => Base::V,
        8 => Base::T,
        9 => Base::W,
        10 => Base::Y,
        11 => Base::H,
        12 => Base::K,
        13 => Base::D,
        14 => Base::B,
        15 => Base::N,
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_sequence() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut src: &[u8], expected: &Sequence) -> io::Result<()> {
            let mut actual = Sequence::default();
            get_sequence(&mut src, &mut actual, expected.len())?;
            assert_eq!(&actual, expected);
            Ok(())
        }

        t(&[], &Sequence::default())?;
        t(&[0x12, 0x40], &"ACG".parse()?)?;
        t(&[0x12, 0x48], &"ACGT".parse()?)?;

        Ok(())
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0), Base::Eq);
        assert_eq!(decode_base(1), Base::A);
        assert_eq!(decode_base(2), Base::C);
        assert_eq!(decode_base(3), Base::M);
        assert_eq!(decode_base(4), Base::G);
        assert_eq!(decode_base(5), Base::R);
        assert_eq!(decode_base(6), Base::S);
        assert_eq!(decode_base(7), Base::V);
        assert_eq!(decode_base(8), Base::T);
        assert_eq!(decode_base(9), Base::W);
        assert_eq!(decode_base(10), Base::Y);
        assert_eq!(decode_base(11), Base::H);
        assert_eq!(decode_base(12), Base::K);
        assert_eq!(decode_base(13), Base::D);
        assert_eq!(decode_base(14), Base::B);
        assert_eq!(decode_base(15), Base::N);
    }
}
