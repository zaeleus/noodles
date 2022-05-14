//! FASTA record sequence base complement.

use std::{error, fmt, iter::FusedIterator, slice};

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
        write!(f, "invalid base: {:04x?}", self.0)
    }
}

impl<'a> Iterator for Complement<'a> {
    type Item = Result<u8, ComplementError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().copied().map(complement)
    }
}

impl<'a> DoubleEndedIterator for Complement<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back().copied().map(complement)
    }
}

impl<'a> ExactSizeIterator for Complement<'a> {}

impl<'a> FusedIterator for Complement<'a> {}

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
        assert_eq!(complement(b'X'), Err(ComplementError(b'X')));
    }
}
