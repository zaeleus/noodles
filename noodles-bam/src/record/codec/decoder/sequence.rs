use std::{error, fmt, num};

use noodles_sam::alignment::record_buf::Sequence;

/// An error when a raw BAM record sequence fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The length is invalid.
    InvalidLength(num::TryFromIntError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
            Self::InvalidLength(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
        }
    }
}

pub(super) fn read_length(src: &mut &[u8]) -> Result<usize, DecodeError> {
    read_u32_le(src).and_then(|n| usize::try_from(n).map_err(DecodeError::InvalidLength))
}

pub(super) fn read_sequence(
    src: &mut &[u8],
    sequence: &mut Sequence,
    base_count: usize,
) -> Result<(), DecodeError> {
    let len = base_count.div_ceil(2);

    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or(DecodeError::UnexpectedEof)?;

    *src = rest;

    let bases = buf
        .iter()
        .flat_map(|&b| [decode_base(b >> 4), decode_base(b)]);

    let dst = sequence.as_mut();
    dst.clear();
    dst.extend(bases);
    dst.truncate(base_count);

    Ok(())
}

fn decode_base(n: u8) -> u8 {
    match n & 0x0f {
        0 => b'=',
        1 => b'A',
        2 => b'C',
        3 => b'M',
        4 => b'G',
        5 => b'R',
        6 => b'S',
        7 => b'V',
        8 => b'T',
        9 => b'W',
        10 => b'Y',
        11 => b'H',
        12 => b'K',
        13 => b'D',
        14 => b'B',
        15 => b'N',
        _ => unreachable!(),
    }
}

fn read_u32_le(src: &mut &[u8]) -> Result<u32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(u32::from_le_bytes(*buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_length() {
        let mut src = &8u32.to_le_bytes()[..];
        assert_eq!(read_length(&mut src), Ok(8));

        let mut src = &[][..];
        assert_eq!(read_length(&mut src), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_sequence() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut src: &[u8], buf: &mut Sequence, expected: &Sequence) -> Result<(), DecodeError> {
            buf.as_mut().clear();
            read_sequence(&mut src, buf, expected.len())?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut sequence = Sequence::default();

        t(&[], &mut sequence, &Sequence::default())?;
        t(&[0x12, 0x40], &mut sequence, &Sequence::from(b"ACG"))?;
        t(&[0x12, 0x48], &mut sequence, &Sequence::from(b"ACGT"))?;

        sequence.as_mut().clear();
        let mut src = &b""[..];
        assert_eq!(
            read_sequence(&mut src, &mut sequence, 4),
            Err(DecodeError::UnexpectedEof)
        );

        Ok(())
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0), b'=');
        assert_eq!(decode_base(1), b'A');
        assert_eq!(decode_base(2), b'C');
        assert_eq!(decode_base(3), b'M');
        assert_eq!(decode_base(4), b'G');
        assert_eq!(decode_base(5), b'R');
        assert_eq!(decode_base(6), b'S');
        assert_eq!(decode_base(7), b'V');
        assert_eq!(decode_base(8), b'T');
        assert_eq!(decode_base(9), b'W');
        assert_eq!(decode_base(10), b'Y');
        assert_eq!(decode_base(11), b'H');
        assert_eq!(decode_base(12), b'K');
        assert_eq!(decode_base(13), b'D');
        assert_eq!(decode_base(14), b'B');
        assert_eq!(decode_base(15), b'N');
    }
}
