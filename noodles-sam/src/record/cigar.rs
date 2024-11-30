use std::{fmt, io, iter};

use crate::{alignment::record::cigar::Op, io::reader::record_buf::cigar::op};

/// Raw SAM record CIGAR operations.
#[derive(Eq, PartialEq)]
pub struct Cigar<'a>(&'a [u8]);

impl<'a> Cigar<'a> {
    /// Creates SAM record CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Cigar;
    /// let cigar = Cigar::new(b"8M13N");
    /// ```
    pub fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::Cigar;
    ///
    /// let cigar = Cigar::new(b"");
    /// assert!(cigar.is_empty());
    ///
    /// let cigar = Cigar::new(b"8M13N");
    /// assert!(!cigar.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns an iterator over CIGAR operations.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     alignment::record::cigar::{op::Kind, Op},
    ///     record::Cigar,
    /// };
    ///
    /// let cigar = Cigar::new(b"");
    /// assert!(cigar.iter().next().is_none());
    ///
    /// let cigar = Cigar::new(b"8M13N");
    /// let mut iter = cigar.iter();
    /// assert_eq!(iter.next().transpose()?, Some(Op::new(Kind::Match, 8)));
    /// assert_eq!(iter.next().transpose()?, Some(Op::new(Kind::Skip, 13)));
    /// assert!(iter.next().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
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

impl fmt::Debug for Cigar<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl crate::alignment::record::Cigar for Cigar<'_> {
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

impl AsRef<[u8]> for Cigar<'_> {
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
