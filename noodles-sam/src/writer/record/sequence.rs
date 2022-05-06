use std::io::{self, Write};

use super::MISSING;
use crate::record::Sequence;

pub fn write_sequence<W>(writer: &mut W, sequence: &Sequence) -> io::Result<()>
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_sequence() -> io::Result<()> {
        use crate::record::{sequence::Base, Sequence};

        let mut buf = Vec::new();
        write_sequence(&mut buf, &Sequence::default())?;
        assert_eq!(buf, b"*");

        buf.clear();
        let sequence = Sequence::from(vec![Base::A, Base::C, Base::G, Base::T]);
        write_sequence(&mut buf, &sequence)?;
        assert_eq!(buf, b"ACGT");

        Ok(())
    }
}
