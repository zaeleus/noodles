use std::{io, mem};

use crate::record::{sequence::Base, Sequence};

pub(crate) fn parse_sequence(src: &[u8], sequence: &mut Sequence) -> io::Result<()> {
    let mut bases = Vec::from(mem::take(sequence));

    for &n in src {
        let base = Base::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        bases.push(base);
    }

    *sequence = Sequence::from(bases);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let mut sequence = Sequence::default();

        sequence.clear();
        parse_sequence(b"", &mut sequence)?;
        let expected = Sequence::default();
        assert_eq!(sequence, expected);

        sequence.clear();
        parse_sequence(b"ACGT", &mut sequence)?;
        let expected = Sequence::from(vec![Base::A, Base::C, Base::G, Base::T]);
        assert_eq!(sequence, expected);

        sequence.clear();
        assert!(matches!(
            parse_sequence(&[0x07], &mut sequence),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
