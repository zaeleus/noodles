mod cigar;
mod data;

pub use self::{cigar::write_cigar, data::write_data};

use std::io::{self, Write};

use crate::record::{QualityScores, Sequence};

const MISSING: u8 = b'*';

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
        for &score in quality_scores.as_ref() {
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
