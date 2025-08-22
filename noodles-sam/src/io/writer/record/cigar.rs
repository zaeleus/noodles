mod op;

use std::io::{self, Write};

use self::op::write_op;
use crate::alignment::record::Cigar;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::{
        record::cigar::{Op, op::Kind},
        record_buf::Cigar as CigarBuf,
    };

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
}
