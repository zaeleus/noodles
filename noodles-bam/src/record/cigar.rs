//! BAM CIGAR and operations.

mod op;

pub use self::op::Op;

use std::{convert::TryFrom, fmt, mem, ops::Deref};

use noodles_sam::record::cigar::op::Kind;

/// BAM record CIGAR.
pub struct Cigar<'a>(&'a [u8]);

impl<'a> Cigar<'a> {
    /// Creates a CIGAR by wrapping raw CIGAR data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Cigar;
    ///
    /// // 36M8S
    /// let data = [0x40, 0x02, 0x00, 0x00, 0x84, 0x00, 0x00, 0x00];
    /// let cigar = Cigar::new(&data);
    ///
    /// assert_eq!(*cigar, data);
    /// ```
    pub fn new(bytes: &[u8]) -> Cigar<'_> {
        Cigar(bytes)
    }

    /// Returns a iterator over the operations in the CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::{cigar::Op, Cigar};
    /// use noodles_sam::record::cigar::op::Kind;
    ///
    /// // 36M8S
    /// let data = [0x40, 0x02, 0x00, 0x00, 0x84, 0x00, 0x00, 0x00];
    /// let cigar = Cigar::new(&data);
    ///
    /// let mut ops = cigar.ops();
    ///
    /// assert_eq!(ops.next(), Some(Op::new(Kind::Match, 36)));
    /// assert_eq!(ops.next(), Some(Op::new(Kind::SoftClip, 8)));
    /// assert_eq!(ops.next(), None);
    /// ```
    pub fn ops(&self) -> Ops<'_> {
        Ops {
            cigar: self.0,
            i: 0,
        }
    }

    /// Calculates the alignment span over the reference sequence.
    ///
    /// This sums the lengths of the CIGAR operations that consume the reference sequence, i.e.,
    /// alignment matches (`M`), deletions from the reference (`D`), skipped reference regions
    /// (`S`), sequence matches (`=`), and sequence mismatches (`X`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::{cigar::Op, Cigar};
    /// use noodles_sam::record::cigar::op::Kind;
    ///
    /// // 36M4D8S
    /// let data = [0x40, 0x02, 0x00, 0x00, 0x43, 0x00, 0x00, 0x00, 0x84, 0x00, 0x00, 0x00];
    /// let cigar = Cigar::new(&data);
    ///
    /// assert_eq!(cigar.mapped_len(), 40);
    /// ```
    pub fn mapped_len(&self) -> u32 {
        self.ops()
            .filter_map(|op| match op.kind() {
                Kind::Match | Kind::Deletion | Kind::Skip | Kind::SeqMatch | Kind::SeqMismatch => {
                    Some(op.len())
                }
                _ => None,
            })
            .sum()
    }
}

impl<'a> fmt::Debug for Cigar<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_list().entries(self.ops()).finish()
    }
}

impl<'a> fmt::Display for Cigar<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for op in self.ops() {
            write!(f, "{}", op)?;
        }

        Ok(())
    }
}

impl<'a> Deref for Cigar<'a> {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        self.0
    }
}

/// An iterator over the operations of a CIGAR.
///
/// This is created by calling [`Cigar::ops`].
///
/// [`Cigar::ops`]: struct.Cigar.html#method.ops
pub struct Ops<'a> {
    cigar: &'a [u8],
    i: usize,
}

impl<'a> Iterator for Ops<'a> {
    type Item = Op;

    fn next(&mut self) -> Option<Self::Item> {
        let size = mem::size_of::<u32>();
        let start = self.i * size;

        if start < self.cigar.len() {
            let end = start + size;

            let data = &self.cigar[start..end];
            let op = Op::try_from(data).unwrap();

            self.i += 1;

            Some(op)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_bytes() {
        let bytes = [0x40, 0x02, 0x00, 0x00, 0x62, 0x03, 0x00, 0x00];
        let cigar = Cigar::new(&bytes);
        let mut ops = cigar.ops();
        assert_eq!(ops.next(), Some(Op::try_from(0x240).unwrap()));
        assert_eq!(ops.next(), Some(Op::try_from(0x362).unwrap()));
        assert_eq!(ops.next(), None);
    }
}
