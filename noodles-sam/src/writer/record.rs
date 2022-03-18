mod data;

pub use self::data::write_data;

use std::io::{self, Write};

use crate::record::{Cigar, QualityScores, Sequence};

const MISSING: u8 = b'*';

pub(super) fn write_cigar<W>(writer: &mut W, cigar: &Cigar) -> io::Result<()>
where
    W: Write,
{
    use crate::record::cigar::op::Kind;

    if cigar.is_empty() {
        writer.write_all(&[MISSING])?;
    } else {
        for op in cigar.iter() {
            write!(writer, "{}", op.len())?;

            let c = match op.kind() {
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

            writer.write_all(&[c])?;
        }
    }

    Ok(())
}

pub(super) fn write_sequence<W>(writer: &mut W, sequence: &Sequence) -> io::Result<()>
where
    W: Write,
{
    if sequence.is_empty() {
        writer.write_all(&[MISSING])?;
    } else {
        for &base in sequence.as_ref() {
            let n = u8::from(base);
            writer.write_all(&[n])?;
        }
    }

    Ok(())
}

pub(super) fn write_quality_scores<W>(
    writer: &mut W,
    quality_scores: &QualityScores,
) -> io::Result<()>
where
    W: Write,
{
    const MIN_VALUE: u8 = b'!';

    if quality_scores.is_empty() {
        writer.write_all(&[MISSING])?;
    } else {
        for &score in quality_scores.iter() {
            let n = u8::from(score) + MIN_VALUE;
            writer.write_all(&[n])?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_cigar() -> io::Result<()> {
        use crate::record::{
            cigar::{op::Kind, Op},
            Cigar,
        };

        let mut buf = Vec::new();
        write_cigar(&mut buf, &Cigar::default())?;
        assert_eq!(buf, b"*");

        buf.clear();
        let cigar = Cigar::from(vec![Op::new(Kind::Match, 8)]);
        write_cigar(&mut buf, &cigar)?;
        assert_eq!(buf, b"8M");

        Ok(())
    }

    #[test]
    fn test_write_sequence() -> io::Result<()> {
        use crate::record::sequence::Base;

        let mut buf = Vec::new();
        write_sequence(&mut buf, &Sequence::default())?;
        assert_eq!(buf, b"*");

        buf.clear();
        let sequence = Sequence::from(vec![Base::A, Base::C, Base::G, Base::T]);
        write_sequence(&mut buf, &sequence)?;
        assert_eq!(buf, b"ACGT");

        Ok(())
    }

    #[test]
    fn test_write_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        write_quality_scores(&mut buf, &QualityScores::default())?;
        assert_eq!(buf, b"*");

        buf.clear();
        let quality_scores = "NDLS".parse()?;
        write_quality_scores(&mut buf, &quality_scores)?;
        assert_eq!(buf, b"NDLS");

        Ok(())
    }
}
