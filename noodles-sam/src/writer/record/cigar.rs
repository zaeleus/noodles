use std::io::{self, Write};

use crate::{
    alignment::record::{cigar::op::Kind, Cigar},
    writer::num,
};

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

            num::write_usize(writer, op.len())?;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Cigar as CigarBuf;

    #[test]
    fn test_write_cigar() -> io::Result<()> {
        use crate::alignment::record::cigar::Op;

        let mut buf = Vec::new();
        write_cigar(&mut buf, &CigarBuf::default())?;
        assert_eq!(buf, b"*");

        buf.clear();
        let cigar: CigarBuf = [Op::new(Kind::Match, 8)].into_iter().collect();
        write_cigar(&mut buf, &cigar)?;
        assert_eq!(buf, b"8M");

        Ok(())
    }
}
