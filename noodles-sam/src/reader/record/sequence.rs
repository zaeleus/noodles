use std::io;

use crate::record::{sequence::Base, Sequence};

pub(crate) fn parse_sequence(src: &[u8]) -> io::Result<Sequence> {
    const MISSING: &[u8] = b"*";

    if src == MISSING {
        return Ok(Sequence::default());
    }

    src.iter()
        .copied()
        .map(Base::try_from)
        .collect::<Result<Vec<_>, _>>()
        .map(Sequence::from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let actual = parse_sequence(b"")?;
        let expected = Sequence::default();
        assert_eq!(actual, expected);

        let actual = parse_sequence(b"ACGT")?;
        let expected = Sequence::from(vec![Base::A, Base::C, Base::G, Base::T]);
        assert_eq!(actual, expected);

        assert!(matches!(
            parse_sequence(&[0x07]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
