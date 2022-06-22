use std::io;

use crate::record::{
    cigar::{op::Kind, Op},
    Cigar,
};

pub(crate) fn parse_cigar(mut src: &[u8]) -> io::Result<Cigar> {
    const MISSING: &[u8] = b"*";

    if src == MISSING {
        return Ok(Cigar::default());
    }

    let mut cigar = Cigar::default();

    while !src.is_empty() {
        let op = parse_op(&mut src)?;
        cigar.as_mut().push(op);
    }

    Ok(cigar)
}

fn parse_op(src: &mut &[u8]) -> io::Result<Op> {
    let len = parse_len(src)?;
    let kind = parse_kind(src)?;
    Ok(Op::new(kind, len))
}

fn parse_len(src: &mut &[u8]) -> io::Result<usize> {
    let (len, i) = lexical_core::parse_partial(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *src = &src[i..];

    Ok(len)
}

fn parse_kind(src: &mut &[u8]) -> io::Result<Kind> {
    let (n, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    match n {
        b'M' => Ok(Kind::Match),
        b'I' => Ok(Kind::Insertion),
        b'D' => Ok(Kind::Deletion),
        b'N' => Ok(Kind::Skip),
        b'S' => Ok(Kind::SoftClip),
        b'H' => Ok(Kind::HardClip),
        b'P' => Ok(Kind::Pad),
        b'=' => Ok(Kind::SequenceMatch),
        b'X' => Ok(Kind::SequenceMismatch),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CIGAR op kind",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar() -> Result<(), Box<dyn std::error::Error>> {
        let src = b"1M13N144S";
        let actual = parse_cigar(src)?;

        let expected = Cigar::try_from(vec![
            Op::new(Kind::Match, 1),
            Op::new(Kind::Skip, 13),
            Op::new(Kind::SoftClip, 144),
        ])?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_len() {
        assert!(matches!(
            parse_len(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
    }

    #[test]
    fn test_parse_kind() -> io::Result<()> {
        fn t(mut src: &[u8], expected: Kind) -> io::Result<()> {
            assert_eq!(parse_kind(&mut src)?, expected);
            Ok(())
        }

        t(b"M", Kind::Match)?;
        t(b"I", Kind::Insertion)?;
        t(b"D", Kind::Deletion)?;
        t(b"N", Kind::Skip)?;
        t(b"S", Kind::SoftClip)?;
        t(b"H", Kind::HardClip)?;
        t(b"P", Kind::Pad)?;
        t(b"=", Kind::SequenceMatch)?;
        t(b"X", Kind::SequenceMismatch)?;

        assert!(matches!(
            parse_kind(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof,
        ));

        assert!(matches!(
            parse_kind(&mut &b"!"[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
