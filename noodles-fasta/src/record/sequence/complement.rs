//! FASTA record sequence base complement.

use std::{error, fmt, iter::FusedIterator, slice};

use bstr::ByteSlice;

/// An iterator that returns the complement of a sequence.
pub struct Complement<'a> {
    iter: slice::Iter<'a, u8>,
}

impl<'a> Complement<'a> {
    pub(super) fn new(iter: slice::Iter<'a, u8>) -> Self {
        Self { iter }
    }
}

/// An error returned when a base does not have a complement.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ComplementError(u8);

impl error::Error for ComplementError {}

impl fmt::Display for ComplementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bytes = [self.0];
        write!(f, "invalid base: {:?}", bytes.as_bstr())
    }
}

impl Iterator for Complement<'_> {
    type Item = Result<u8, ComplementError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().copied().map(complement)
    }
}

impl DoubleEndedIterator for Complement<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back().copied().map(complement)
    }
}

impl ExactSizeIterator for Complement<'_> {}

impl FusedIterator for Complement<'_> {}

fn complement(b: u8) -> Result<u8, ComplementError> {
    match b {
        b'A' => Ok(b'T'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'T' => Ok(b'A'),
        b'U' => Ok(b'A'),
        b'W' => Ok(b'W'),
        b'S' => Ok(b'S'),
        b'M' => Ok(b'K'),
        b'K' => Ok(b'M'),
        b'R' => Ok(b'Y'),
        b'Y' => Ok(b'R'),
        b'B' => Ok(b'V'),
        b'D' => Ok(b'H'),
        b'H' => Ok(b'D'),
        b'V' => Ok(b'B'),
        b'N' => Ok(b'N'),

        b'a' => Ok(b't'),
        b'c' => Ok(b'g'),
        b'g' => Ok(b'c'),
        b't' => Ok(b'a'),
        b'u' => Ok(b'a'),
        b'w' => Ok(b'w'),
        b's' => Ok(b's'),
        b'm' => Ok(b'k'),
        b'k' => Ok(b'm'),
        b'r' => Ok(b'y'),
        b'y' => Ok(b'r'),
        b'b' => Ok(b'v'),
        b'd' => Ok(b'h'),
        b'h' => Ok(b'd'),
        b'v' => Ok(b'b'),
        b'n' => Ok(b'n'),

        _ => Err(ComplementError(b)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() -> Result<(), ComplementError> {
        let complement = Complement::new(b"ACGT".iter());
        let actual: Vec<_> = complement.collect::<Result<_, _>>()?;
        let expected = b"TGCA";
        assert_eq!(actual, expected);
        Ok(())
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), Ok(b'T'));
        assert_eq!(complement(b'C'), Ok(b'G'));
        assert_eq!(complement(b'G'), Ok(b'C'));
        assert_eq!(complement(b'T'), Ok(b'A'));
        assert_eq!(complement(b'U'), Ok(b'A'));
        assert_eq!(complement(b'W'), Ok(b'W'));
        assert_eq!(complement(b'S'), Ok(b'S'));
        assert_eq!(complement(b'M'), Ok(b'K'));
        assert_eq!(complement(b'K'), Ok(b'M'));
        assert_eq!(complement(b'R'), Ok(b'Y'));
        assert_eq!(complement(b'Y'), Ok(b'R'));
        assert_eq!(complement(b'B'), Ok(b'V'));
        assert_eq!(complement(b'D'), Ok(b'H'));
        assert_eq!(complement(b'H'), Ok(b'D'));
        assert_eq!(complement(b'V'), Ok(b'B'));
        assert_eq!(complement(b'N'), Ok(b'N'));

        assert_eq!(complement(b'a'), Ok(b't'));
        assert_eq!(complement(b'c'), Ok(b'g'));
        assert_eq!(complement(b'g'), Ok(b'c'));
        assert_eq!(complement(b't'), Ok(b'a'));
        assert_eq!(complement(b'u'), Ok(b'a'));
        assert_eq!(complement(b'w'), Ok(b'w'));
        assert_eq!(complement(b's'), Ok(b's'));
        assert_eq!(complement(b'm'), Ok(b'k'));
        assert_eq!(complement(b'k'), Ok(b'm'));
        assert_eq!(complement(b'r'), Ok(b'y'));
        assert_eq!(complement(b'y'), Ok(b'r'));
        assert_eq!(complement(b'b'), Ok(b'v'));
        assert_eq!(complement(b'd'), Ok(b'h'));
        assert_eq!(complement(b'h'), Ok(b'd'));
        assert_eq!(complement(b'v'), Ok(b'b'));
        assert_eq!(complement(b'n'), Ok(b'n'));

        assert_eq!(complement(b'X'), Err(ComplementError(b'X')));
    }
}
