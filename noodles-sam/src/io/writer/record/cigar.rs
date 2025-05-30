use std::io::{self, Write};

use crate::{
    alignment::record::{
        Cigar,
        cigar::{Op, op::Kind},
    },
    io::writer::num,
};

/// Writes a SAM record CIGAR string.
///
/// # Examples
///
/// ```
/// use noodles_sam::{
///     alignment::{
///         record::cigar::{op::Kind, Op},
///         record_buf::Cigar,
///     },
///     io::writer::record::write_cigar,
/// };
///
/// let mut buf = Vec::new();
/// let cigar = Cigar::default();
/// write_cigar(&mut buf, &cigar)?;
/// assert_eq!(buf, b"*");
///
/// let mut buf = Vec::new();
/// let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
/// write_cigar(&mut buf, &cigar)?;
/// assert_eq!(buf, b"4M");
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn write_cigar<W, C>(writer: &mut W, cigar: &C) -> io::Result<()>
where
    W: Write,
    C: Cigar,
{
    use super::MISSING;

    if cigar.is_empty() {
        writer.write_all(&[MISSING])?;
    } else {
        for result in cigar.iter() {
            let op = result?;
            write_op(writer, op)?;
        }
    }

    Ok(())
}

fn write_op<W>(writer: &mut W, op: Op) -> io::Result<()>
where
    W: Write,
{
    num::write_usize(writer, op.len())?;
    write_kind(writer, op.kind())?;
    Ok(())
}

fn write_kind<W>(writer: &mut W, kind: Kind) -> io::Result<()>
where
    W: Write,
{
    let c = match kind {
        Kind::Match => b'M',
        Kind::Insertion => b'I',
        Kind::Deletion => b'D',
        Kind::Skip => b'N',
        Kind::SoftClip => b'S',
        Kind::HardClip => b'H',
        Kind::Pad => b'P',
        Kind::SequenceMatch => b'=',
        Kind::SequenceMismatch => b'X',
    };

    writer.write_all(&[c])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::record_buf::Cigar as CigarBuf;

    #[test]
    fn test_write_cigar() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, cigar: &CigarBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_cigar(buf, cigar)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let cigar = CigarBuf::default();
        t(&mut buf, &cigar, b"*")?;

        let cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
        t(&mut buf, &cigar, b"4M")?;

        let cigar: CigarBuf = [Op::new(Kind::Match, 4), Op::new(Kind::HardClip, 2)]
            .into_iter()
            .collect();
        t(&mut buf, &cigar, b"4M2H")?;

        Ok(())
    }

    #[test]
    fn test_write_op() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_op(&mut buf, Op::new(Kind::Match, 1))?;
        assert_eq!(buf, b"1M");

        buf.clear();
        write_op(&mut buf, Op::new(Kind::Match, 1 << 28))?;
        assert_eq!(buf, b"268435456M");

        Ok(())
    }

    #[test]
    fn test_write_kind() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, kind: Kind, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_kind(buf, kind)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Kind::Match, b"M")?;
        t(&mut buf, Kind::Insertion, b"I")?;
        t(&mut buf, Kind::Deletion, b"D")?;
        t(&mut buf, Kind::Skip, b"N")?;
        t(&mut buf, Kind::SoftClip, b"S")?;
        t(&mut buf, Kind::HardClip, b"H")?;
        t(&mut buf, Kind::Pad, b"P")?;
        t(&mut buf, Kind::SequenceMatch, b"=")?;
        t(&mut buf, Kind::SequenceMismatch, b"X")?;

        Ok(())
    }
}
