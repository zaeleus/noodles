use std::{io, iter};

use crate::{alignment::record::cigar::Op, io::reader::record_buf::cigar::op};

/// Raw SAM record CIGAR operations.
#[derive(Debug, Eq, PartialEq)]
pub struct Cigar<'a>(&'a [u8]);

impl<'a> Cigar<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any CIGAR operations.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns an iterator over CIGAR operations.
    pub fn iter(&self) -> impl Iterator<Item = Result<Op, op::ParseError>> + '_ {
        use crate::io::reader::record_buf::cigar::op::parse_op;

        let mut src = self.0;

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(parse_op(&mut src))
            }
        })
    }
}

impl<'a> crate::alignment::record::field::Cigar for Cigar<'a> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.as_ref()
            .iter()
            .filter(|&b| {
                matches!(
                    b,
                    b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'=' | b'X'
                )
            })
            .count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Op>> + '_> {
        Box::new(self.iter().map(|result| {
            result
                .map(|op| Op::new(op.kind(), op.len()))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        }))
    }
}

impl<'a> AsRef<[u8]> for Cigar<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Cigar<'a>> for crate::alignment::record_buf::Cigar {
    type Error = io::Error;

    fn try_from(Cigar(src): Cigar<'a>) -> Result<Self, Self::Error> {
        use crate::io::reader::record_buf::parse_cigar;

        let mut cigar = Self::default();

        if !src.is_empty() {
            parse_cigar(src, &mut cigar)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        }

        Ok(cigar)
    }
}

#[cfg(test)]
mod tests {
    use crate::alignment::record::cigar::op::Kind;

    use super::*;

    #[test]
    fn test_iter() -> Result<(), op::ParseError> {
        let cigar = Cigar::new(b"");
        assert!(cigar.iter().next().is_none());

        let cigar = Cigar::new(b"8M13N");
        let actual: Vec<_> = cigar.iter().collect::<Result<_, _>>()?;
        let expected = [Op::new(Kind::Match, 8), Op::new(Kind::Skip, 13)];
        assert_eq!(actual, expected);

        Ok(())
    }
}
